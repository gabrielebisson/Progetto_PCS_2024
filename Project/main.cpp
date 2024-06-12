#include "Utils.hpp"
#include "readFileTest.hpp"
#include <iostream>
#include <fstream>
#include "readFileTest.hpp"
#include "StructDFN.hpp"
#include "Eigen/Eigen"

using namespace LibraryDFN;
using namespace Eigen;

int main()
{
    double tol=__DBL_EPSILON__;
    DFN dfn;
    dfn.numFratture=1;
    dfn.idFratture={0};
    dfn.numVertici={{7}};
    dfn.vertici={{Vector3d(5.,2,0.)
                   ,Vector3d(11.,0.,0.)
                   ,Vector3d(20.,0.,0.)
                   ,Vector3d(22.,8.,0.)
                   ,Vector3d(14.,16.,0.)
                   ,Vector3d(2.,13.,0.)
                    ,Vector3d(0.,7.,0.)}};
    dfn.versori={Vector3d(0.,0.,1.)};

    dfn.idTracce={0,1,2,3,4,5,6};
    dfn.estremiTracce={{Vector3d(1.,10.,0.),Vector3d(17.,13.,0.)}
                      ,{Vector3d(6.,14.,0.),Vector3d(8.,1.,0.)}
                      ,{Vector3d(19.,11.,0.),Vector3d(19.,0.,0.)}
                      ,{Vector3d(10.,13.,0.),Vector3d(6.,6.,0.)}
                      ,{Vector3d(15.,9.,0.),Vector3d(19.,11.,0.)}
                      ,{Vector3d(15.,6.,0.),Vector3d(19.,5.,0.)}
                      ,{Vector3d(11.,7.,0.),Vector3d(12.,5.,0.)}};

    dfn.traccePassanti={{0,1,2}};
    dfn.tracceNonPassanti={{3,4,5,6}};
    definisci_mesh(dfn,tol);
    // PolygonalMesh mesh;
    // mesh.NumberCell2D=5;
    // mesh.Cell2DVertices={{0,4,8,9,5,1},{2,5,9,7},{7,10,11,6,3},{9,11,10},{9,8,6}};
    // mesh.Cell0DCoordinates= {Vector3d(2.,1.,0.),Vector3d(2.,-1.5,0.),Vector3d(-1.5,-3.5,0.),Vector3d(-4.5,1.,0.),Vector3d(1.,3.,0.),Vector3d(0.66,-2.27,0.),Vector3d(-2.48,1.74,0.),Vector3d(-2.41,-2.13,0.),Vector3d(-0.27,2.54,0.),Vector3d(-1.32,0.25,0.),Vector3d(-1.87,-0.96,0.),Vector3d(-1.94,1.05,0.)};
    // Vector3d est1(-1.5,-3.5,0.);
    // Vector3d est2(-1.32,0.25,0.);
    // unsigned int id=0;
    // std::tuple<std::vector<unsigned int>,std::vector<std::array<unsigned int,2>>> zio=contatto_poligoni_segmento(dfn,id,mesh,est1,est2,tol);
    return 0;

}
