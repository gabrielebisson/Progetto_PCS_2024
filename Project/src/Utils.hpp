#pragma once

#include "StructDFN.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include <list>
#include <tuple>
#include <cmath>

namespace LibraryDFN
{
    //funzioni relative all'input/output
    void readDFNFromFile(const std::string& filename, DFN& dfn);
    void printTraces(const DFN& dfn, const std::string& filename);
    void printTracesByFracture(const DFN& dfn, const std::string& filename);

    //funzioni comuni al calcolo delle tracce e alla creazione delle mesh poligonali
    std::vector<std::array<unsigned int,3>> triangola_frattura(DFN& disc_frac_net,unsigned int idr);
    inline Eigen::Matrix<double,3,3> allinea_xy(const std::vector<Eigen::Vector3d>& vertici, Eigen::Vector3d normale);

    //funzioni relative al calcolo delle tracce
    inline Eigen::Vector3d versore_normale(const std::vector<Eigen::Vector3d>& poligono);
    std::vector<std::array<unsigned int,2>> scarta_fratture(DFN& disc_frac_net,const std::vector<Eigen::Vector3d>& versori_normali,const double& tol);
    void memorizza_tracce(DFN& disc_frac_net,double tol);

    //funzioni relative alla creazione delle mesh poligonali
    std::tuple<Eigen::Vector3d,bool> interseca_segmenti(const Eigen::Vector3d& A1,const Eigen::Vector3d& A2,const Eigen::Vector3d& B1,const Eigen::Vector3d& B2,const double& tol);
    inline Eigen::Vector3d versore_normale_2(const DFN& disc_frac_net,const std::vector<unsigned int>& poligono,const unsigned int& pos);
    std::array<bool,2> interno_bordo_poligono(const DFN& disc_frac_net, const std::vector<unsigned int>& poligono,const Eigen::Vector3d& x,const unsigned int& pos_frattura, const double& tol);
}
