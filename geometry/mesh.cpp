#include "mesh.h"
#include "element.h"
#include "cohesivezone.h"
#include "czinsertiontool.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <ios>
#include <iomanip>
#include <iostream>
#include <map>

#include <Eigen/Core>


icy::ConcurrentPool<icy::Node> icy::Mesh::NodeFactory(reserveConst);
icy::ConcurrentPool<icy::Element> icy::Mesh::ElementFactory(reserveConst);
icy::ConcurrentPool<icy::CohesiveZone> icy::Mesh::CZFactory(reserveConst);


icy::Mesh::Mesh()
{
}

icy::Mesh::~Mesh()
{
    Reset();
}

void icy::Mesh::Reset()
{
    NodeFactory.release(nodes);
    ElementFactory.release(elems);
    CZFactory.release(czs);
}



icy::Node* icy::Mesh::AddNode()
{
    Node* nd = NodeFactory.take();
    nd->Reset();
    nd->globId = (int)nodes.size();
    nodes.push_back(nd);
    return nd;
}

icy::Element* icy::Mesh::AddElement()
{
    Element* elem = ElementFactory.take();
    elem->Reset();
    elem->elemId = (int)elems.size();
    elems.push_back(elem);
    return elem;
}

icy::CohesiveZone* icy::Mesh::AddCZ()
{
    CohesiveZone *cz = CZFactory.take();
    cz->Reset();
    czs.push_back(cz);
    return cz;
}


void icy::Mesh::LoadMSH(const std::string &fileName)
{
    Reset();
    gmsh::clear();
    gmsh::open(fileName);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    std::unordered_map<std::size_t, std::size_t> mtags; // gmsh nodeTag -> sequential position in nodes[]

    // GET NODES
    gmsh::model::mesh::getNodesByElementType(4, nodeTags, nodeCoords, parametricCoords);


    NodeFactory.release(nodes);
    ElementFactory.release(elems);
    nodes.clear();
    elems.clear();

    // set the size of the resulting nodes array
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in deformable");

        Node *nd = AddNode();
        mtags[tag] = nd->globId;
        nd->x0 = nd->xn = Eigen::Vector3d(nodeCoords[i*3+0], nodeCoords[i*3+1], nodeCoords[i*3+2]);
    }

    // GET ELEMENTS - per grain (entity)
    std::vector<std::pair<int,int>> dimTagsGrains;
    gmsh::model::getEntities(dimTagsGrains,3);

    for(std::size_t j=0;j<dimTagsGrains.size();j++)
    {
        std::vector<std::size_t> tetraTags, nodeTagsInTetra;
        int entityTag = dimTagsGrains[j].second;
        gmsh::model::mesh::getElementsByType(4, tetraTags, nodeTagsInTetra,entityTag);

        for(std::size_t i=0;i<tetraTags.size();i++)
        {
            icy::Element *elem = AddElement();
            elem->grainId = (int)j;
            for(int k=0;k<4;k++) elem->nds[k] = nodes[mtags.at(nodeTagsInTetra[i*4+k])];
        }
    }

    std::vector<std::size_t> lineTags, nodeTagsInLines;
    gmsh::model::mesh::getElementsByType(1, lineTags, nodeTagsInLines);
    std::cout << "number of lines: " << lineTags.size() << std::endl;


    // center mesh
    Eigen::Vector3d offset(0.5,0.5,0);

    for(Node *nd : nodes)
    {
        nd->x0 -= offset;
        nd->xn = nd->xt = nd->x0;
        nd->pinned = nd->x0.z() < 1e-7;
    }


    CZInsertionTool czit;
    czit.InsertCZs(*this);

    for(icy::Element *elem : elems) elem->Precompute();     // Dm matrix and volume
    MarkIncidentFaces();

    gmsh::clear();
}


void icy::Mesh::MarkIncidentFaces()
{
    // initialize incident faces information in nodes and elements

    // (1) find exposed faces
    std::map<std::tuple<int,int,int>,icy::CZInsertionTool::Facet> facets;
    for(icy::Element *e : elems)
    {
        for(int k=0;k<4;k++)
        {
            std::tuple<int,int,int> key = icy::CZInsertionTool::Facet::make_key(
                        e->nds[Element::fi[k][0]],
                    e->nds[Element::fi[k][1]],
                    e->nds[Element::fi[k][2]]);
            icy::CZInsertionTool::Facet facet;
            facet.key = key;
            facet.elems[0] = e;
            facet.facet_idx[0] = k;
            auto result = facets.insert({key,facet});
            if(result.second == false)
            {
                icy::CZInsertionTool::Facet &f = result.first->second;
                f.elems[1] = e;
                f.facet_idx[1] = k;
            }
        }
    }

    // (2) distribute exposed faces to their incident nodes
    for(Node *nd : nodes) nd->incident_faces.clear();

    int count = 0;
    for(auto &kvp : facets)
    {
        icy::CZInsertionTool::Facet &f = kvp.second;
        if(f.elems[1] != nullptr) continue;
        Element *e = f.elems[0];
        unsigned k = f.facet_idx[0];

        uint32_t facet_code = (uint32_t)f.facet_idx[0] | (uint32_t)f.elems[0]->elemId << 2;
        for(int i=0;i<3;i++)
        {
            Node *nd = e->nds[Element::fi[k][i]];
            nd->incident_faces.push_back(facet_code);
        }
        count++;
    }

    // (3) per element - collect all exposed faces from element's 4 nodes
    for(Element *elem : elems)
    {
        auto &r = elem->elem_incident_faces;
        r.clear();
        for(int i=0;i<4;i++)
            for(uint32_t facet_code : elem->nds[i]->incident_faces)
                r.push_back(facet_code);
        std::sort(r.begin(),r.end());
        r.resize(std::distance(r.begin(),std::unique(r.begin(), r.end())));
    }
}











void icy::Mesh::ExportForAbaqus(std::string fileName, double czStrength, std::string jobName, std::string batchName,
                                double YoungsModulus, double czElasticity, double czEnergy)
{
    std::ofstream s;
    s.open(fileName,std::ios_base::trunc|std::ios_base::out);
    s << std::setprecision(9);
    s << "from abaqus import *\n";
    s << "from abaqusConstants import *\n";
    s << "from caeModules import *\n";

    s << "import mesh\n";
    s << "import regionToolset\n";
    s << "import os\n";

    QString path = QDir::currentPath()+ "/" + QString::fromStdString(batchName) + "/inp";

    s << "os.chdir(r\"" << path.toStdString() << "\")\n";

    s << "p = mdb.models['Model-1'].Part(name='MyPart1', dimensionality=THREE_D, type=DEFORMABLE_BODY)\n";

    for(Node *nd : nodes)
        s << "p.Node(coordinates=(" << nd->xn[0] << "," << nd->xn[1] << "," << nd->xn[2] << "))\n";

    s << "n = p.nodes\n";

    for(Element *e : elems)
        s << "p.Element(nodes=(n["<<e->nds[0]->globId<<"],n["<<e->nds[1]->globId<<
             "],n["<<e->nds[3]->globId<<"],n["<<e->nds[2]->globId<<"]), elemShape=TET4)\n";

    for(icy::CohesiveZone *c : czs)
        s << "p.Element(nodes=(n["<<c->nds[0]->globId<<
             "], n["<<c->nds[1]->globId<<
             "], n["<<c->nds[2]->globId<<
             "], n["<<c->nds[3]->globId<<
             "], n["<<c->nds[4]->globId<<
             "], n["<<c->nds[5]->globId<<"]), elemShape=WEDGE6)\n";

    s << "elemType_bulk = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD, secondOrderAccuracy=OFF, distortionControl=DEFAULT)\n";
    s << "elemType_coh = mesh.ElemType(elemCode=COH3D6, elemLibrary=STANDARD)\n";

    // region1 - bulk elements

    s << "region1 = p.elements[0:"<<elems.size()<<"]\n";
    s << "p.setElementType(regions=(region1,), elemTypes=(elemType_bulk,))\n";
    s << "p.Set(elements=(region1,), name='Set-1-elems')\n";

    s << "region2cz = p.elements[" << elems.size() << ":" << elems.size()+czs.size() << "]\n";
    s << "p.setElementType(regions=(region2cz,), elemTypes=(elemType_coh,))\n";
    s << "p.Set(elements=(region2cz,), name='Set-2-czs')\n";

    // region2 - pinned nodes
    s << "region3pinned = (";
    for(Node *nd : nodes)
        if(nd->pinned)
            s << "p.nodes["<<nd->globId<<":"<<nd->globId+1<<"],";
    s << ")\n";
    s << "p.Set(nodes=region3pinned,name='Set3-pinned')\n";

    // create bulk material
    s << "mat1 = mdb.models['Model-1'].Material(name='Material-1-bulk')\n";
    s << "mat1.Density(table=((900.0, ), ))\n";
    s << "mat1.Elastic(table=((" << YoungsModulus << ", 0.3), ))\n";

    // cz material
    s << "mat2 = mdb.models['Model-1'].Material(name='Material-2-czs')\n";
    s << "mat2.Density(table=((1.0, ), ))\n";
    s << "mat2.MaxsDamageInitiation(table=((" << czStrength << "," << czStrength/2 << "," << czStrength << "), ))\n";
    s << "mat2.maxsDamageInitiation.DamageEvolution(type=ENERGY, table=((" << czEnergy << ", ), ))\n";
    s << "mat2.Elastic(type=TRACTION, table=((" << czElasticity << "," << czElasticity/2 << "," << czElasticity/2 << "), ))\n";

    // sections
    s << "mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1-bulk', "
        "material='Material-1-bulk', thickness=None)\n";
    s << "mdb.models['Model-1'].CohesiveSection(name='Section-2-czs', "
        "material='Material-2-czs', response=TRACTION_SEPARATION, "
        "outOfPlaneThickness=None)\n";

    // section assignments
    s << "region = p.sets['Set-1-elems']\n";
    s << "p.SectionAssignment(region=region, sectionName='Section-1-bulk', offset=0.0, "
        "offsetType=MIDDLE_SURFACE, offsetField='', "
        "thicknessAssignment=FROM_SECTION)\n";
    s << "region = p.sets['Set-2-czs']\n";
    s << "p = mdb.models['Model-1'].parts['MyPart1']\n";
    s << "p.SectionAssignment(region=region, sectionName='Section-2-czs', offset=0.0, "
        "offsetType=MIDDLE_SURFACE, offsetField='', "
        "thicknessAssignment=FROM_SECTION)\n";

    // indenter
    s << "s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2.0)\n";
    s << "g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n";
    s << "s.setPrimaryObject(option=STANDALONE)\n";
    s << "s.ArcByCenterEnds(center=(0.0, 0.0), point1=(-0.05, 0.0), point2=(0.05, -0.0125), direction=COUNTERCLOCKWISE)\n";
    s << "p2 = mdb.models['Model-1'].Part(name='Part-2', dimensionality=THREE_D, type=ANALYTIC_RIGID_SURFACE)\n";
    //s << "p2 = mdb.models['Model-1'].parts['Part-2']\n";
    s << "p2.AnalyticRigidSurfExtrude(sketch=s, depth=1.0)\n";
    s << "s.unsetPrimaryObject()\n";

    s << "v1 = p2.vertices\n";
    s << "p2.ReferencePoint(point=v1[2])\n";


    // assembly
    s << "a1 = mdb.models['Model-1'].rootAssembly\n";
    s << "a1.DatumCsysByDefault(CARTESIAN)\n";

    // add and rotate main part
    s << "inst1 = a1.Instance(name='MyPart1-1', part=p, dependent=ON)\n";
    s << "a1.rotate(instanceList=('MyPart1-1', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0), angle=-90.0)\n";

    // add and rotate indenter
    s << "a1.Instance(name='Part-2-1', part=p2, dependent=ON)\n";
    s << "a1.translate(instanceList=('Part-2-1', ), vector=(0.0, 1.05, 0.0))\n";


    // create step
    s << "mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=2.0, improvedDtMethod=ON)\n";

    // create field output request
    s << "mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(numIntervals=300)\n";

    // gravity load
    s << "mdb.models['Model-1'].Gravity(name='Load-1', createStepName='Step-1',comp2=-10.0, distributionType=UNIFORM, field='')\n";

    // BC - pinned nodes
    s << "region = inst1.sets['Set3-pinned']\n";
    s << "mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Initial', region=region, localCsys=None)\n";

    // BC - moving indenter
    s << "r1 = a1.instances['Part-2-1'].referencePoints\n";
    s << "refPoints1=(r1[2], )\n";
    s << "region = a1.Set(referencePoints=refPoints1, name='Set-1-indenterRP')\n";
    s << "mdb.models['Model-1'].VelocityBC(name='BC-2', createStepName='Step-1', "
        "region=region, v1=0.0, v2=-0.001, v3=0.0, vr1=0.0, vr2=0.0, vr3=0.0, "
        "amplitude=UNSET, localCsys=None, distributionType=UNIFORM, fieldName='')\n";

    // rigid body constraint
    s << "s1 = a1.instances['Part-2-1'].faces\n";
    s << "side2Faces1 = s1[0:1]\n";
    s << "region5=a1.Surface(side2Faces=side2Faces1, name='Surf-1')\n";
    s << "r1 = a1.instances['Part-2-1'].referencePoints\n";
    s << "refPoints1=(r1[2], )\n";
    s << "region1=regionToolset.Region(referencePoints=refPoints1)\n";
    s << "mdb.models['Model-1'].RigidBody(name='Constraint-1', refPointRegion=region1, surfaceRegion=region5)\n";

    // create interaction property
    s << "mdb.models['Model-1'].ContactProperty('IntProp-1')\n";
    s << "mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior("
        "formulation=FRICTIONLESS)\n";
    s << "mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior("
        "pressureOverclosure=HARD, allowSeparation=ON, "
        "constraintEnforcementMethod=DEFAULT)\n";

    // create interaction itself
    s << "mdb.models['Model-1'].ContactExp(name='Int-1', createStepName='Step-1')\n";
    s << "mdb.models['Model-1'].interactions['Int-1'].includedPairs.setValuesInStep("
        "stepName='Step-1', useAllstar=ON)\n";
    s << "mdb.models['Model-1'].interactions['Int-1'].contactPropertyAssignments.appendInStep("
        "stepName='Step-1', assignments=((GLOBAL, SELF, 'IntProp-1'), ))\n";


    // record indenter force

    s << "mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-1', name="
        "'H-Output-2', rebar=EXCLUDE, region="
        "mdb.models['Model-1'].rootAssembly.sets['Set-1-indenterRP'], sectionPoints="
        "DEFAULT, timeInterval=0.0001, variables=('RF2', ))\n";

    //create job
    s << "mdb.Job(name='" << jobName << "', model='Model-1', description='', type=ANALYSIS,"
    "atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,"
    "memoryUnits=PERCENTAGE, explicitPrecision=SINGLE,"
    "nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF,"
    "contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='',"
    "resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=4,"
    "activateLoadBalancing=False, numThreadsPerMpiProcess=1,"
    "multiprocessingMode=DEFAULT, numCpus=4)\n";

    // write .inp file
    s << "mdb.jobs['" << jobName << "'].writeInput(consistencyChecking=OFF)";

    s.close();
}

void icy::Mesh::RotateSample(double angleInDegrees)
{
    if(angleInDegrees == 0) return;
    double alpha = angleInDegrees*M_PI/180.;
    double cosA = std::cos(alpha);
    double sinA = std::sin(alpha);
    Eigen::Matrix3d R;
    R << cosA, -sinA, 0,
            sinA, cosA, 0,
            0, 0, 1;
    for(Node *nd : nodes)
    {
        nd->xn = nd->xt = R*nd->x0;
    }
}
