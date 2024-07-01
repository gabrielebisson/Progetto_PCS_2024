#pragma once

#include <vector>
#include <map>
#include <array>
#include <list>
#include "Eigen/Eigen"

namespace PolygonalLibrary {
    struct PolygonalMesh
    {
        //celle0D
        unsigned int NumberCell0D = 0; //numero celle 0D (punti)
        std::vector<unsigned int> Cell0DId = {}; //vettore con gli ID dei punti
        std::vector<Eigen::Vector3d> Cell0DCoordinates = {}; //vettore con le coordinate dei punti (vettore di array con x e y)
        std::vector<unsigned int> Cell0DMarkers = {}; //vettore con i marker dei vari punti

        //celle1D
        unsigned int NumberCell1D = 0; //numero celle 1D (segmenti)
        std::vector<unsigned int> Cell1DId = {}; //vettore con gli ID dei segmenti
        std::vector<std::array<unsigned int,2>> Cell1DVertices = {}; //vettore con gli ID dei punti che fanno da estremi al segmento (ID_inizio,ID_fine)
        std::vector<unsigned int> Cell1DMarkers = {}; //vettore con i marker dei vari lati

        //celle2D
        unsigned int NumberCell2D = 0; //numero celle 2D (poligoni)
        std::vector<unsigned int> Cell2DId = {}; //vettore con gli ID dei poligoni
        std::vector<std::vector<unsigned int>> Cell2DVertices = {}; //vettore di vettori con gli ID dei vertici dei poligoni (in senso antiorario)
        std::vector<std::vector<unsigned int>> Cell2DEdges = {}; //vettore di vettori con gli ID dei lati dei poligoni (in senso antiorario)
    };
}
