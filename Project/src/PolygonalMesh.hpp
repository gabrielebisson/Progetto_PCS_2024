#pragma once

#include <vector>
#include <map>
#include <array>
#include <list>

using namespace std;

namespace PolygonalLibrary {
    struct PolygonalMesh
    {
        //celle0D
        unsigned int NumberCell0D = 0; //numero celle 0D (punti)
        vector<unsigned int> Cell0DId = {}; //vettore con gli ID dei punti
        vector<array<double,2>> Cell0DCoordinates = {}; //vettore con le coordinate dei punti (vettore di array con x e y)
        map<unsigned int,list<unsigned int>> Cell0DMarkers = {}; //dizionario con chiave il marker e con argomento la lista dei punti che hanno quel marker

        //celle1D
        unsigned int NumberCell1D = 0; //numero celle 1D (segmenti)
        vector<unsigned int> Cell1DId = {}; //vettore con gli ID dei segmenti
        vector<array<unsigned int,2>> Cell1DVertices = {}; //vettore con gli ID dei punti che fanno da estremi al segmento (ID_inizio,ID_fine)
        map<unsigned int,list<unsigned int>> Cell1DMarkers = {}; //dizionario con chiave il marker e con argomento la lista dei segmenti che hanno quel marker

        //celle2D
        unsigned int NumberCell2D = 0; //numero celle 2D (poligoni)
        vector<unsigned int> Cell2DId = {}; //vettore con gli ID dei poligoni
        vector<vector<unsigned int>> Cell2DVertices = {}; //vettore di vettori con gli ID dei vertici dei poligoni (in senso antiorario)
        vector<vector<unsigned int>> Cell2DEdges = {}; //vettore di vettori con gli ID dei lati dei poligoni (in senso antiorario)
    };

    //Questa struttura serve solamente come supporto per la mesh poligonale
    struct MeshDiSupporto
    {
        //celle0D
        unsigned int NumberCell0D = 0;
        list<unsigned int> Cell0DId = {};
        list<array<double,2>> Cell0DCoordinates = {};
        map<unsigned int,list<unsigned int>> Cell0DMarkers = {};

        //celle1D
        unsigned int NumberCell1D = 0;
        list<unsigned int> Cell1DId = {};
        list<array<unsigned int,2>> Cell1DVertices = {};
        map<unsigned int,list<unsigned int>> Cell1DMarkers = {};

        //celle2D
        unsigned int NumberCell2D = 0;
        list<unsigned int> Cell2DId = {};
        list<list<unsigned int>> Cell2DVertices = {};
        list<list<unsigned int>> Cell2DEdges = {};
    };
}

