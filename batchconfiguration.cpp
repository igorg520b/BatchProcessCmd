#include "batchconfiguration.h"

#include <QFileInfo>
#include <QTextStream>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <map>

#include <spdlog/spdlog.h>

#include "mesh2d.h"
#include "mesh.h"
#include "element.h"
#include "element2d.h"
#include "cohesivezone.h"
#include "cohesivezone2d.h"


QString BatchConfiguration::BatchName() const
{
    return ToShortName(batchFileName);
}


void BatchConfiguration::Load(QString fileName)
{
    QFile f(fileName);
    f.open(QIODevice::ReadOnly | QIODevice::Text);
    if(!f.exists()) qDebug() << "config file does not exist";
    QTextStream ts(&f);

    if (!f.exists())
    {
        qDebug() << "config file not open";
        throw std::runtime_error("");
    }

    QMap<QString, QString> kvp;
    QString line;
    while(!ts.atEnd())
    {
        line = ts.readLine();
        QStringList split = line.split("=");
        kvp.insert(split[0],split[1]);
    }
    f.close();

    if(kvp.contains("dim")) dim = kvp.value("dim").toInt();



    mshFileName = kvp.value("mesh");
    if(kvp.contains("czsStrength")) czsStrength = kvp.value("czStrength").toDouble();
    if(kvp.contains("YoungsModulus")) YoungsModulus = kvp.value("YoungsModulus").toDouble();

    if(kvp.contains("czElasticity")) czElasticity = kvp.value("czElasticity").toDouble();
    if(kvp.contains("czEnergy")) czEnergy = kvp.value("czEnergy").toDouble();

    numberOfCores = kvp.value("nCPUs").toInt();

    if(kvp.contains("indenterRadius")) indenterRadius = kvp.value("indenterRadius").toDouble();
    if(kvp.contains("indenterDepth")) indenterDepth = kvp.value("indenterDepth").toDouble();
    if(kvp.contains("indentationRate")) indentationRate = kvp.value("indentationRate").toDouble();
    if(kvp.contains("horizontalOffset")) horizontalOffset = kvp.value("horizontalOffset").toDouble();

    if(kvp.contains("insertCZs")) insertCZs = kvp.value("insertCZs")!="false";

    if(kvp.contains("confinement"))
    {
        QString s = kvp.value("confinement");
        if(s=="bottomOnly") confinement=ConfinementType::bottomOnly;
        else if(s=="frontAndBack") confinement=ConfinementType::frontAndBack;
        else if(s=="full") confinement=ConfinementType::full;
        else if(s=="sides") confinement=ConfinementType::sides;
    }
//{1: ^{0}d}|{2: ^{0}.3e}
    spdlog::info("test");
    spdlog::info("{0:<15} : {1}", "mshFileName", mshFileName.toStdString().c_str());
    spdlog::info("{0:<15} : {1}", "dim", dim);
    spdlog::info("{0:<15} : {1:03.2e}", "czsStrength", czsStrength);
    spdlog::info("{0:<15} : {1:03.2e}", "YoungsModulus", YoungsModulus);
    spdlog::info("{0:<15} : {1:03.2e}", "czElasticity", czElasticity);
    spdlog::info("{0:<15} : {1}", "czEnergy", czEnergy);
    spdlog::info("{0:<15} : {1}", "numberOfCores", numberOfCores);
    spdlog::info("{0:<15} : {1}", "indenterRadius", indenterRadius);
    spdlog::info("{0:<15} : {1}", "indenterDepth", indenterDepth);
    spdlog::info("{0:<15} : {1}", "indentationRate", indentationRate);
    spdlog::info("{0:<15} : {1}", "horizontalOffset", horizontalOffset);
    spdlog::info("{0:<15} : {1}", "insertCZs", insertCZs);
    spdlog::info("{0:<15} : {1}", "confinement", confinement);

    //
    batchFileName = fileName;
}



QString BatchConfiguration::ToShortName(QString filePath)
{
    QFileInfo fi(filePath);
    return fi.baseName();
}

void BatchConfiguration::GeneratePythonScript()
{
    spdlog::info("BatchConfiguration::GeneratePythonScript()");

    // create directory for .py file
    QString tentativeDir = QDir::currentPath()+ "/" + BatchName();
    QDir dir(tentativeDir);
    if(!dir.exists())
    {
        spdlog::info("creating directory {}",tentativeDir.toStdString().c_str());
        dir.mkdir(tentativeDir);
    }

    if(dim==2)
    {
        // 2D MESH
        qDebug() << "BatchConfiguration::ProducePYFiles - 2D";
        icy::Mesh2D m2d;
        QString meshPath = "meshes2d/" + mshFileName;
        QFileInfo f(meshPath);
        if(!f.exists())
        {
            qDebug() << "mesh file not found" << meshPath;
            throw std::runtime_error("mesh file does not exist");
        }

        QString pyPath = QDir::currentPath()+ "/" + BatchName() + "/" + BatchName()+".py";
        qDebug() << "loading 2d mesh " << meshPath;
        m2d.LoadMSH(meshPath.toStdString(), insertCZs);


        qDebug() << "ExportForAbaqus - 2D";

        std::ofstream s;
        s.open(pyPath.toStdString(), std::ios_base::trunc|std::ios_base::out);
        s << std::setprecision(9);
        s << "from abaqus import *\n";
        s << "from abaqusConstants import *\n";
        s << "from caeModules import *\n";

        s << "import mesh\n";
        s << "import regionToolset\n";
        s << "import os\n";

        s << "p = mdb.models['Model-1'].Part(name='MyPart1', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)\n";

        for(icy::Node2D *nd : m2d.nodes)
            s << "p.Node(coordinates=(" << nd->x0[0] << "," << nd->x0[1] << ",0))\n";

        s << "n = p.nodes\n";

        for(icy::Element2D *e : m2d.elems)
            s << "p.Element(nodes=(n["<<e->nds[0]->globId<<"],n["<<e->nds[1]->globId<<
                 "],n["<<e->nds[2]->globId<<"]), elemShape=TRI3)\n";



        if(m2d.czs.size()>0)
            for(icy::CohesiveZone2D *c : m2d.czs)
                s << "p.Element(nodes=(n["<<c->nds[0]->globId<<
                     "], n["<<c->nds[1]->globId<<
                     "], n["<<c->nds[3]->globId<<
                     "], n["<<c->nds[2]->globId<<
                     "]), elemShape=QUAD4)\n";

        s << "elemType_bulk = mesh.ElemType(elemCode=CPS3, elemLibrary=STANDARD, secondOrderAccuracy=OFF, distortionControl=DEFAULT)\n";

        if(m2d.czs.size()>0)
            s << "elemType_coh = mesh.ElemType(elemCode=COH2D4, elemLibrary=STANDARD)\n";


        // region1 - bulk elements
        s << "region1 = p.elements[0:" << m2d.elems.size() << "]\n";
        s << "p.setElementType(regions=(region1,), elemTypes=(elemType_bulk,))\n";
        s << "p.Set(elements=(region1,), name='Set-1-elems')\n";

        if(m2d.czs.size()>0)
        {
            s << "region2cz = p.elements[" << m2d.elems.size() << ":" << m2d.elems.size() + m2d.czs.size() << "]\n";
            s << "p.setElementType(regions=(region2cz,), elemTypes=(elemType_coh,))\n";
            s << "p.Set(elements=(region2cz,), name='Set-2-czs')\n";
        }

        // region - pinned nodes

        for(icy::Node2D *nd : m2d.nodes)
            if(nd->x0.y() < 1e-7) nd->pinned = true;


        if(confinement == 1)
        {
            // full
            for(icy::Node2D *nd : m2d.nodes)
                if(nd->x0.y() < 0.5 && (nd->x0.x() < 1e-7|| nd->x0.x() > 2.5-1e-7)) nd->pinned = true;
        }



        s << "region3pinned = (";
        for(icy::Node2D *nd : m2d.nodes)
            if(nd->pinned)
                s << "p.nodes["<<nd->globId<<":"<<nd->globId+1<<"],";
        s << ")\n";

        s << "p.Set(nodes=region3pinned,name='Set3-pinned')\n";


        // create bulk material
        s << "mat1 = mdb.models['Model-1'].Material(name='Material-1-bulk')\n";
        s << "mat1.Density(table=((900.0, ), ))\n";
        s << "mat1.Elastic(table=((" << YoungsModulus << ", 0.3), ))\n";

        // cz material
        if(m2d.czs.size()>0)
        {
            s << "mat2 = mdb.models['Model-1'].Material(name='Material-2-czs')\n";
            s << "mat2.Density(table=((1.0, ), ))\n";
            s << "mat2.MaxsDamageInitiation(table=((" << czsStrength << "," << czsStrength/2 << "," << czsStrength/2 << "), ))\n";
            s << "mat2.maxsDamageInitiation.DamageEvolution(type=ENERGY, table=((" << czEnergy << ", ), ))\n";
            s << "mat2.Elastic(type=TRACTION, table=((" << czElasticity << "," << czElasticity/2 << "," << czElasticity/2 << "), ))\n";
        }

        // sections
        s << "mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1-bulk', "
             "material='Material-1-bulk', thickness=None)\n";

        if(m2d.czs.size()>0)
            s << "mdb.models['Model-1'].CohesiveSection(name='Section-2-czs', "
                 "material='Material-2-czs', response=TRACTION_SEPARATION, "
                 "outOfPlaneThickness=None)\n";


        // section assignments
        s << "region = p.sets['Set-1-elems']\n";
        s << "p.SectionAssignment(region=region, sectionName='Section-1-bulk', offset=0.0, "
             "offsetType=MIDDLE_SURFACE, offsetField='', "
             "thicknessAssignment=FROM_SECTION)\n";

        if(m2d.czs.size()>0)
        {
            s << "region = p.sets['Set-2-czs']\n";
            s << "p = mdb.models['Model-1'].parts['MyPart1']\n";
            s << "p.SectionAssignment(region=region, sectionName='Section-2-czs', offset=0.0, "
                 "offsetType=MIDDLE_SURFACE, offsetField='', "
                 "thicknessAssignment=FROM_SECTION)\n";
        }


        // indenter
        s << "s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2.0)\n";
        s << "g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n";
        s << "s.ArcByCenterEnds(center=(0.0, 0.0), point1=("<< -indenterRadius << ", 0.0), point2=(" << indenterRadius << ", -0.0125), direction=COUNTERCLOCKWISE)\n";

        s << "p2 = mdb.models['Model-1'].Part(name='Part-2', dimensionality=TWO_D_PLANAR, type=ANALYTIC_RIGID_SURFACE)\n";

        s << "p2.AnalyticRigidSurf2DPlanar(sketch=s)\n";


        s << "v1 = p2.vertices\n";

        s << "p2.ReferencePoint(point=p2.InterestingPoint(p2.edges[0], CENTER))\n";


        // assembly
        s << "a1 = mdb.models['Model-1'].rootAssembly\n";
        s << "a1.DatumCsysByDefault(CARTESIAN)\n";

        // add and rotate main part
        s << "inst1 = a1.Instance(name='MyPart1-1', part=p, dependent=ON)\n";

        horizontalOffset += -sqrt(pow(indenterRadius,2)-pow(indenterRadius-indenterDepth,2))-1e-7;

        double xOffset = horizontalOffset;
        double yOffset = -indenterDepth + indenterRadius + 1;
        double zOffset = 0;
        s << "a1.Instance(name='Part-2-1', part=p2, dependent=ON)\n";
        // rotate indenter
        s << "a1.rotate(instanceList=('Part-2-1', ), axisPoint=(0.0, 0.0, 0.0)," << "axisDirection=(0.0, 0.0, 1.0), angle=45.0)\n";
        s << "a1.translate(instanceList=('Part-2-1', ), vector=("<< xOffset << ", " << yOffset << ", " << zOffset << "))\n";


        // create step
        double timePeriod = 10;
        int numIntervals = 200*timePeriod;
        s << "mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=" << timePeriod << ", improvedDtMethod=ON)\n";

        // create field output request
        s << "mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(numIntervals=" << numIntervals << ")\n";

        // gravity load
        s << "mdb.models['Model-1'].Gravity(name='Load-1', createStepName='Step-1',comp2=-10.0, distributionType=UNIFORM, field='')\n";

        // BC - pinned nodes
        s << "region = inst1.sets['Set3-pinned']\n";
        s << "mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Initial', region=region, localCsys=None)\n";

        // BC - moving indenter
        double xVelocity = indentationRate;
        double yVelocity = 0;
        s << "r1 = a1.instances['Part-2-1'].referencePoints\n";
        s << "refPoints1=(r1[2], )\n";
        s << "region = a1.Set(referencePoints=refPoints1, name='Set-1-indenterRP')\n";
        s << "mdb.models['Model-1'].VelocityBC(name='BC-2', createStepName='Step-1', "
             "region=region, v1="<< xVelocity << ", v2=" << yVelocity << ", v3=0.0, vr1=0.0, vr2=0.0, vr3=0.0, "
                                                                         "amplitude=UNSET, localCsys=None, distributionType=UNIFORM, fieldName='')\n";



        s << "ed1 = a1.instances['Part-2-1'].edges\n";
        s << "side2Edges1 = ed1[0:1]\n";
        s << "a1.Surface(name='Surf-1', side2Edges=side2Edges1)\n";

        s << "mdb.models['Model-1'].RigidBody(name='Constraint-1', refPointRegion=regionToolset.Region(referencePoints="
             "(a1.instances['Part-2-1'].referencePoints[2],)), surfaceRegion=a1.surfaces['Surf-1'])\n";


        // rigid body constraint

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
             "DEFAULT, timeInterval=0.0001, variables=('RF1','RF2', ))\n";

        //create job
        s << "mdb.Job(name='" << BatchName().toStdString() << "', model='Model-1', description='', type=ANALYSIS,"
                                            "atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,"
                                            "memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE,"
                                            "nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,"
                                            "contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='',"
                                            "resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains="
          <<numberOfCores<<","
                           "activateLoadBalancing=False, numThreadsPerMpiProcess=1,"
                           "multiprocessingMode=DEFAULT, numCpus="<<numberOfCores<<")\n";

        s.close();
        qDebug() << "ExportForAbaqus done";


    }
    else if(dim==3)
    {
        // 3D MESH
        qDebug() << "BatchConfiguration::ProducePYFiles - 3D";
        icy::Mesh m;

        QString meshPath = "meshes/" + mshFileName;
        QFileInfo f(meshPath);
        if(!f.exists())
        {
            qDebug() << "mesh file not found" << meshPath;
            throw std::runtime_error("mesh file does not exist");
        }

        QString pyPath = QDir::currentPath()+ "/" + BatchName() + "/" + BatchName()+".py";
        qDebug() << "loading 3d mesh " << meshPath;
        m.LoadMSH(meshPath.toStdString(), insertCZs);

        qDebug() << "ExportForAbaqus";

        std::ofstream s;
        s.open(pyPath.toStdString(), std::ios_base::trunc|std::ios_base::out);
        s << std::setprecision(9);
        s << "from abaqus import *\n";
        s << "from abaqusConstants import *\n";
        s << "from caeModules import *\n";

        s << "import mesh\n";
        s << "import regionToolset\n";
        s << "import os\n";

        s << "p = mdb.models['Model-1'].Part(name='MyPart1', dimensionality=THREE_D, type=DEFORMABLE_BODY)\n";

        for(icy::Node *nd : m.nodes)
            s << "p.Node(coordinates=(" << nd->x0[0] << "," << nd->x0[1] << "," << nd->x0[2] << "))\n";

        s << "n = p.nodes\n";

        for(icy::Element *e : m.elems)
            s << "p.Element(nodes=(n["<<e->nds[0]->globId<<"],n["<<e->nds[1]->globId<<
                 "],n["<<e->nds[3]->globId<<"],n["<<e->nds[2]->globId<<"]), elemShape=TET4)\n";


        if(m.czs.size()>0)
            for(icy::CohesiveZone *c : m.czs)
                s << "p.Element(nodes=(n["<<c->nds[0]->globId<<
                     "], n["<<c->nds[1]->globId<<
                     "], n["<<c->nds[2]->globId<<
                     "], n["<<c->nds[3]->globId<<
                     "], n["<<c->nds[4]->globId<<
                     "], n["<<c->nds[5]->globId<<"]), elemShape=WEDGE6)\n";

        s << "elemType_bulk = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD, secondOrderAccuracy=OFF, distortionControl=DEFAULT)\n";

        if(m.czs.size()>0)
            s << "elemType_coh = mesh.ElemType(elemCode=COH3D6, elemLibrary=STANDARD)\n";

        // region1 - bulk elements
        s << "region1 = p.elements[0:" << m.elems.size() << "]\n";
        s << "p.setElementType(regions=(region1,), elemTypes=(elemType_bulk,))\n";
        s << "p.Set(elements=(region1,), name='Set-1-elems')\n";

        if(m.czs.size()>0)
        {
            s << "region2cz = p.elements[" << m.elems.size() << ":" << m.elems.size() + m.czs.size() << "]\n";
            s << "p.setElementType(regions=(region2cz,), elemTypes=(elemType_coh,))\n";
            s << "p.Set(elements=(region2cz,), name='Set-2-czs')\n";
        }

        // region - pinned nodes
        if(confinement == 1)
        {
            // full
            for(icy::Node *nd : m.nodes)
                if(nd->x0.z() < 0.5 && (nd->x0.x() < 1e-7 || nd->x0.y() < 1e-7 ||
                                        nd->x0.x() > 2.5-1e-7 || nd->x0.y() > 1.5-1e-5)) nd->pinned = true;
        }
        else if(confinement == 3)
        {
            // front and back
            for(icy::Node *nd : m.nodes)
                if(nd->x0.z() < 0.5 && (nd->x0.x() < 1e-7 ||
                                        nd->x0.x() > 2.5-1e-7)) nd->pinned = true;
        }
        else if(confinement == 2)
        {
            // sides
            for(icy::Node *nd : m.nodes)
                if(nd->x0.z() < 0.5 && (nd->x0.y() < 1e-7 ||
                                        nd->x0.y() > 1.5-1e-5)) nd->pinned = true;
        }

        s << "region3pinned = (";
        for(icy::Node *nd : m.nodes)
            if(nd->pinned)
                s << "p.nodes["<<nd->globId<<":"<<nd->globId+1<<"],";
        s << ")\n";

        s << "p.Set(nodes=region3pinned,name='Set3-pinned')\n";

        // create bulk material
        s << "mat1 = mdb.models['Model-1'].Material(name='Material-1-bulk')\n";
        s << "mat1.Density(table=((900.0, ), ))\n";
        s << "mat1.Elastic(table=((" << YoungsModulus << ", 0.3), ))\n";

        // cz material
        if(m.czs.size()>0)
        {
            s << "mat2 = mdb.models['Model-1'].Material(name='Material-2-czs')\n";
            s << "mat2.Density(table=((1.0, ), ))\n";
            s << "mat2.MaxsDamageInitiation(table=((" << czsStrength << "," << czsStrength/2 << "," << czsStrength/2 << "), ))\n";
            s << "mat2.maxsDamageInitiation.DamageEvolution(type=ENERGY, table=((" << czEnergy << ", ), ))\n";
            s << "mat2.Elastic(type=TRACTION, table=((" << czElasticity << "," << czElasticity/2 << "," << czElasticity/2 << "), ))\n";
        }

        // sections
        s << "mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1-bulk', "
             "material='Material-1-bulk', thickness=None)\n";

        if(m.czs.size()>0)
            s << "mdb.models['Model-1'].CohesiveSection(name='Section-2-czs', "
                 "material='Material-2-czs', response=TRACTION_SEPARATION, "
                 "outOfPlaneThickness=None)\n";

        // section assignments
        s << "region = p.sets['Set-1-elems']\n";
        s << "p.SectionAssignment(region=region, sectionName='Section-1-bulk', offset=0.0, "
             "offsetType=MIDDLE_SURFACE, offsetField='', "
             "thicknessAssignment=FROM_SECTION)\n";

        if(m.czs.size()>0)
        {
            s << "region = p.sets['Set-2-czs']\n";
            s << "p = mdb.models['Model-1'].parts['MyPart1']\n";
            s << "p.SectionAssignment(region=region, sectionName='Section-2-czs', offset=0.0, "
                 "offsetType=MIDDLE_SURFACE, offsetField='', "
                 "thicknessAssignment=FROM_SECTION)\n";
        }

        // indenter
        double indenterLength = 1.5;
        s << "s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2.0)\n";
        s << "g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n";
        s << "s.setPrimaryObject(option=STANDALONE)\n";
        s << "s.ArcByCenterEnds(center=(0.0, 0.0), point1=("<< -indenterRadius << ", 0.0), point2=(" << indenterRadius << ", -0.0125), direction=COUNTERCLOCKWISE)\n";
        s << "p2 = mdb.models['Model-1'].Part(name='Part-2', dimensionality=THREE_D, type=ANALYTIC_RIGID_SURFACE)\n";
        s << "p2.AnalyticRigidSurfExtrude(sketch=s, depth=" << indenterLength << ")\n";
        s << "s.unsetPrimaryObject()\n";

        s << "v1 = p2.vertices\n";
        s << "p2.ReferencePoint(point=v1[2])\n";


        // assembly
        s << "a1 = mdb.models['Model-1'].rootAssembly\n";
        s << "a1.DatumCsysByDefault(CARTESIAN)\n";

        // add and rotate main part
        s << "inst1 = a1.Instance(name='MyPart1-1', part=p, dependent=ON)\n";
        s << "a1.rotate(instanceList=('MyPart1-1', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0), angle=-90.0)\n";

        // add and rotate the indenter
        if(horizontalOffset==0)
        {
            horizontalOffset = -sqrt(pow(indenterRadius,2)-pow(indenterRadius-indenterDepth,2))-1e-7;
        }

        double xOffset = horizontalOffset;
        double yOffset = -indenterDepth + indenterRadius + 1;
        double zOffset = -indenterLength/2;
        s << "a1.Instance(name='Part-2-1', part=p2, dependent=ON)\n";
        // rotate indenter
        s << "a1.rotate(instanceList=('Part-2-1', ), axisPoint=(0.0, 0.0, 0.0),"
             "axisDirection=(0.0, 0.0, 1.0), angle=45.0)\n";
        s << "a1.translate(instanceList=('Part-2-1', ), vector=("<< xOffset << ", " << yOffset << ", " << zOffset << "))\n";


        // create step
        double timePeriod = 10;
        int numIntervals = 200*timePeriod;
        s << "mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=" << timePeriod << ", improvedDtMethod=ON)\n";

        // create field output request
        s << "mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(numIntervals=" << numIntervals << ")\n";

        // gravity load
        s << "mdb.models['Model-1'].Gravity(name='Load-1', createStepName='Step-1',comp2=-10.0, distributionType=UNIFORM, field='')\n";

        // BC - pinned nodes
        s << "region = inst1.sets['Set3-pinned']\n";
        s << "mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Initial', region=region, localCsys=None)\n";

        // BC - moving indenter
        double xVelocity = indentationRate;
        double yVelocity = 0;
        s << "r1 = a1.instances['Part-2-1'].referencePoints\n";
        s << "refPoints1=(r1[2], )\n";
        s << "region = a1.Set(referencePoints=refPoints1, name='Set-1-indenterRP')\n";
        s << "mdb.models['Model-1'].VelocityBC(name='BC-2', createStepName='Step-1', "
             "region=region, v1="<< xVelocity << ", v2=" << yVelocity << ", v3=0.0, vr1=0.0, vr2=0.0, vr3=0.0, "
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
             "DEFAULT, timeInterval=0.0001, variables=('RF1','RF2', ))\n";

        //create job
        s << "mdb.Job(name='" << BatchName().toStdString() << "', model='Model-1', description='', type=ANALYSIS,"
                                            "atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,"
                                            "memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE,"
                                            "nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,"
                                            "contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='',"
                                            "resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains="
          <<numberOfCores<<","
                           "activateLoadBalancing=False, numThreadsPerMpiProcess=1,"
                           "multiprocessingMode=DEFAULT, numCpus="<<numberOfCores<<")\n";

        // write .inp file
//        s << "mdb.jobs['" << BatchName().toStdString() << "'].writeInput(consistencyChecking=OFF)";

        s.close();
        qDebug() << "ExportForAbaqus done";
    }



}

/*

void icy::Mesh::ExportForAbaqus(std::string fileName, double czStrength, std::string jobName, std::string batchName,
                                double YoungsModulus, double czElasticity, double czEnergy,
                                bool rhitaSetup, double indenterRadius, double indenterDepth,
                                double indentationRate, double horizontalOffset, int nCPUs,
                                int confinement)
{


}


*/

