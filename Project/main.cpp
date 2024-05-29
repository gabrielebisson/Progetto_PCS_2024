#include "Utils.hpp"
#include "readFileTest.hpp"
#include <iostream>
#include <fstream>
#include <gtest/gtest.h>
#include "readFileTest.hpp"
#include "StructDFN.hpp"
#include "Eigen/Eigen"

using namespace LibraryDFN;
using namespace Eigen;

int main()
{
    double tol=__DBL_EPSILON__;
    // DFN prova;
    // std::string nome_file="DFN/FR200_data.txt";
    // readDFNFromFile(nome_file,prova);
    // memorizza_tracce(prova,tol);
    // printTraces(prova,"risultato_tracce.txt");
    // printTracesByFracture(prova,"risultato_tracce_ordinate.txt");
    Vector3d A1(-4.,3.,0.);
    Vector3d A2(-2.,5.,4.);
    Vector3d B1(-1.,4.,3.);
    Vector3d B2=B1+2.666*(A1+0.666*(A2-A1)-B1);
    std::tuple<Vector3d,bool> prova=interseca_segmenti(A1,A2,B1,B2,tol);
    return 0;

}
