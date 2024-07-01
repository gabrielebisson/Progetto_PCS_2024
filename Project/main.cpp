#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <gtest/gtest.h>

#include "StructDFN.hpp"
#include "Eigen/Eigen"

using namespace LibraryDFN;
using namespace Eigen;

int main()
{
    double tol=__DBL_EPSILON__;
    DFN prova;
    std::string nome_file="DFN/FR10_data.txt";
    readDFNFromFile(nome_file,prova);
    memorizza_tracce(prova,tol);
    printTraces(prova,"Tracce.txt");
    sortTracesAndPrintByFracture(prova,"Tracce_per_frattura.txt");
    definisci_mesh(prova,tol);
    return 0;
}
