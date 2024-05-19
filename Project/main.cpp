#include "StructDFN.hpp"
#include "Utils.hpp"
#include "Eigen/Eigen"

using namespace LibraryDFN;
using namespace Eigen;

int main()
{
    double tol=__DBL_EPSILON__;
    DFN prova;
    prova.numFratture=8;
    prova.idFratture={0,1,2,3,4,5,6,7};
    prova.numVertici={4,3,3,4,5,3,4,3};
    prova.vertici={{Vector3d(2.5,2.,0.),Vector3d(-3.5,0.5,0.),Vector3d(-5.,4.5,0.),Vector3d(-1.,4.5,0.)},{Vector3d(-4.,3.,0.2),Vector3d(-1.5,4.,2.),Vector3d(0.,0.,-1.)},{Vector3d(0.,0.,-1),Vector3d(-4.,-3.,0.),Vector3d(-1.,-3.5,0.)},{Vector3d(-4.,-3.,0.),Vector3d(-1.,-3.5,0.),Vector3d(-1.,-3.5,4.),Vector3d(-4.,-3.,3.)},{Vector3d(4.5,-2.5,0.),Vector3d(3.5,-5.,0.),Vector3d(5.,-7.,0.),Vector3d(9.,-4.,0.),Vector3d(7.,-1.,0.)},{Vector3d(1.5,-8.5,-3.),Vector3d(3.,-7.5,4.),Vector3d(7.,-3.5,0.)},{Vector3d(5.,0.,0.),Vector3d(7.,0.,0.),Vector3d(7.,2.,0.),Vector3d(5.,2.,0.)},{Vector3d(6.,1.,0.),Vector3d(6.,3.,1.),Vector3d(3.,1.5,7.)}};
    memorizza_tracce(prova,tol);
    return 0;
}
