#pragma once

#include <vector>
#include "Eigen/Eigen"
#include "PolygonalMesh.hpp"

using namespace PolygonalLibrary;

namespace LibraryDFN{

struct DFN {
    //fratture
    unsigned int numFratture; //numero delle fratture
    std::vector<unsigned int> idFratture; //vettore degli id delle fratture
    std::vector<unsigned int> numVertici; //vettore di numeri di vertici di ogni frattura (ordinati nel modo numVertici[i]= numero vertici della frattura idFratture[i]) (se non mi sbaglio serve a qualcosa, in lettura, non in scrittura)
    std::vector<std::vector<Eigen::Vector3d>> vertici; //vettore di matrici di dimensione 3x(numero vertici), ogni colonna individua le coordinate x, y, z del punto (ordinati nel modo vertici[i]= insieme delle coordinate dei vertici della frattura idFratture[i] in senso antiorario)
    std::vector<Eigen::Vector3d> versori; //vettore con tutti i versori normali alle fratture

    //tracce
    unsigned int numTracce; //numero delle tracce
    std::vector<unsigned int> idTracce; //vettore degli id delle tracce
    std::vector<std::array<unsigned int,2>> tracce; //vettore di coppie di fratture che generano la traccia
    std::vector<std::array<Eigen::Vector3d,2>> estremiTracce; //vettore di estremi delle tracce come coordinate 3D
    std::vector<std::array<bool,2>> tips; //vettore di coppie di booleani, ogni booleano indica se la traccia è passante per la frattura che lo genera (vero se non passante, falso se passante)
    std::vector<double> lunghezze; //vettore di lunghezze delle tracce
    std::vector<std::vector<unsigned int>> tracceNonPassanti; //tracce non passanti ordinate per frattura e lunghezza decrescente
    std::vector<std::vector<unsigned int>> traccePassanti; //tracce passanti ordinate per frattura e lunghezza decrescente

    //mesh poligonali
    std::vector<PolygonalMesh> meshPoligonali; //vettore con tutte le mesh poligonali relative ad ogni frattura
};
}
