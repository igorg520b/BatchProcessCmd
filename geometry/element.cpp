#include "element.h"
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include <complex>
#include <cmath>


void icy::Element::Reset()
{
    elem_incident_faces.clear();
    elem_incident_faces.reserve(32);
}


void icy::Element::Precompute()
{
    strain_energy_density = 0;
    principal_stress[0] = 0;
    principal_stress[1] = 0;
    principal_stress[2] = 0;
    max_shear_stress = 0;
    hydrostatic_stress = 0;

    computeDm();
    volume_initial = Dm.determinant()/6.;
    if(volume_initial==0) throw std::runtime_error("element's initial area is zero");
    else if(volume_initial < 0)
    {
        std::swap(nds[0],nds[1]);
        computeDm();
        volume_initial = Dm.determinant()/6.;
    }

    // deformed shape matrix
    Eigen::Matrix3d Ds;
    Ds << (nds[0]->xt - nds[3]->xt), (nds[1]->xt - nds[3]->xt), (nds[2]->xt - nds[3]->xt);
    if(Ds.determinant() <= 0)
    {
        throw std::runtime_error("Ds.determinant");
    }
}


Eigen::Matrix3d icy::Element::getF_at_n() { return getDs_at_n()*DmInv; }

Eigen::Matrix3d icy::Element::getDs_at_n()
{
    return (Eigen::Matrix3d() <<
                (nds[0]->xn - nds[3]->xn),
            (nds[1]->xn - nds[3]->xn),
            (nds[2]->xn - nds[3]->xn)).finished();
}

void icy::Element::computeDm()
{
    Dm << (nds[0]->x0 - nds[3]->x0), (nds[1]->x0 - nds[3]->x0), (nds[2]->x0 - nds[3]->x0);
    DmInv = Dm.inverse();
}


Eigen::Matrix3d icy::Element::getSpecialInverseDs()
{
    // deformed shape matrix
    Eigen::Matrix3d Ds;
    Ds << (nds[1]->xt - nds[0]->xt), (nds[2]->xt - nds[0]->xt), (nds[3]->xt - nds[0]->xt);
    return Ds.inverse();
}

bool icy::Element::IsNodeInside(Eigen::Matrix3d &b, Eigen::Vector3d p)
{
    const double eps = 1e-10;
    Eigen::Vector3d r = b*(p-nds[0]->xt);
    return (r[0]>eps && r[1]>eps && r[2]>eps && r.sum()<(1.-eps));
}




Eigen::Matrix3d icy::Element::DDs[12] = {
    (Eigen::Matrix3d() <<
         1,0,0,
         0,0,0,
         0,0,0).finished(),

    (Eigen::Matrix3d() <<
         0,0,0,
         1,0,0,
         0,0,0).finished(),

    (Eigen::Matrix3d() <<
         0,0,0,
         0,0,0,
         1,0,0).finished(),

    //2
    (Eigen::Matrix3d() <<
         0,1,0,
         0,0,0,
         0,0,0).finished(),

    (Eigen::Matrix3d() <<
         0,0,0,
         0,1,0,
         0,0,0).finished(),

    (Eigen::Matrix3d() <<
         0,0,0,
         0,0,0,
         0,1,0).finished(),

    //3
    (Eigen::Matrix3d() <<
         0,0,1,
         0,0,0,
         0,0,0).finished(),

    (Eigen::Matrix3d() <<
         0,0,0,
         0,0,1,
         0,0,0).finished(),

    (Eigen::Matrix3d() <<
         0,0,0,
         0,0,0,
         0,0,1).finished(),


    //4
    (Eigen::Matrix3d() <<
         -1,-1,-1,
         0,0,0,
         0,0,0).finished(),

    (Eigen::Matrix3d() <<
         0,0,0,
         -1,-1,-1,
         0,0,0).finished(),

    (Eigen::Matrix3d() <<
         0,0,0,
         0,0,0,
         -1,-1,-1).finished()
        };


Eigen::Matrix<double,12,12> icy::Element::consistentMassMatrix =
    (Eigen::Matrix<double,12,12>() <<
         2,0,0, 1,0,0, 1,0,0, 1,0,0,
         0,2,0, 0,1,0, 0,1,0, 0,1,0,
         0,0,2, 0,0,1, 0,0,1, 0,0,1,

         1,0,0, 2,0,0, 1,0,0, 1,0,0,
         0,1,0, 0,2,0, 0,1,0, 0,1,0,
         0,0,1, 0,0,2, 0,0,1, 0,0,1,

         1,0,0, 1,0,0, 2,0,0, 1,0,0,
         0,1,0, 0,1,0, 0,2,0, 0,1,0,
         0,0,1, 0,0,1, 0,0,2, 0,0,1,

         1,0,0, 1,0,0, 1,0,0, 2,0,0,
         0,1,0, 0,1,0, 0,1,0, 0,2,0,
         0,0,1, 0,0,1, 0,0,1, 0,0,2
     ).finished()*(1./20.);


