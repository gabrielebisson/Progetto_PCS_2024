#pragma once

#include "StructDFN.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include <list>
#include <tuple>
#include <cmath>

#include <algorithm>
#include <map>

namespace LibraryDFN
{

    //funzioni relative all'input/output
    void readDFNFromFile(const std::string& filename, DFN& dfn);
    void printTraces(const DFN& dfn, const std::string& filename);
    void sortTracesAndPrintByFracture(DFN &dfn, const std::string& filename);

    //funzioni relative al calcolo delle tracce
    std::vector<std::array<unsigned int,3>> triangola_frattura(DFN& disc_frac_net,unsigned int idr);
    inline Eigen::Matrix<double,3,3> allinea_xy(const std::vector<Eigen::Vector3d>& vertici, Eigen::Vector3d normale);
    inline Eigen::Vector3d versore_normale(const std::vector<Eigen::Vector3d>& poligono);
    std::vector<std::array<unsigned int,2>> scarta_fratture(DFN& disc_frac_net,const double& tol);
    void memorizza_tracce(DFN& disc_frac_net,double tol);

    //funzioni relative alla creazione delle mesh poligonali
    std::tuple<Eigen::Vector3d,double> interseca_segmenti(const Eigen::Vector3d& A1,const Eigen::Vector3d& A2,const Eigen::Vector3d& B1,const Eigen::Vector3d& B2);
    std::tuple<std::vector<unsigned int>,std::vector<std::array<unsigned int,2>>> contatto_poligoni_segmento(const DFN& disc_frac_net,const unsigned int& id,PolygonalMesh& mesh,const Eigen::Vector3d& est1,const Eigen::Vector3d& est2, const double& tol);
    std::array<std::vector<unsigned int>,2> nuovo_poligono(PolygonalMesh& mesh,std::map<unsigned int, std::vector<unsigned int>>& mappa_vecchi_lati_nuovi_lati,std::map<std::array<unsigned int,2>,unsigned int>& mappa_estremi_nuovi_lati_nuovi_lati,const std::vector<unsigned int>& poligoni,const unsigned int& i,const unsigned int& vecchio_lato_partenza,const unsigned int& vecchio_lato_arrivo,const unsigned int& est_taglio_partenza,const unsigned int& est_taglio_arrivo,const unsigned int& ID_1D);
    void definisci_mesh(DFN& disc_frac_net, double tol);
    void aggiorna_mesh(PolygonalMesh& mesh, const std::vector<unsigned int>& poligoni, const std::vector<std::array<unsigned int, 2>> &lati_coinvolti, const Eigen::Vector3d& est1, const Eigen::Vector3d& est2, const double& tol);


}
