#pragma once

#include "StructDFN.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include <list>
#include <tuple>

namespace LibraryDFN
{
    void readDFNFromFile(const std::string& filename, DFN& dfn);
    void printTraces(const DFN& dfn, const std::string& filename);
    void printTracesByFracture(const DFN& dfn, const std::string& filename);
    std::vector<std::array<unsigned int,3>> triangola_frattura(DFN& disc_frac_net,unsigned int idr);
    std::vector<std::array<unsigned int,2>> scarta_fratture(DFN& disc_frac_net);
    Eigen::Matrix<double,3,3> allinea_xy(const std::vector<Eigen::Vector3d>& vertici);
    void memorizza_tracce(DFN& disc_frac_net,double tol);
}