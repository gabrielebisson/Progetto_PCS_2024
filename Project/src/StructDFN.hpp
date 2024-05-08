// definiamo la struttura dati della DFN

#pragma once

#include <vector>

namespace LibraryDFN{
    struct DFN {
        unsigned int numFratture;
        std::vector<unsigned int> idFratture;

        std::vector<unsigned int> numVertici;
        std::vector<std::vector<std::array<double, 3>>> vertici; // x, y, z

        unsigned int numTracce;
        std::vector<unsigned int> idTracce;
        std::vector<std::array<unsigned int,2>> tracce; // idFratture1, idFratture2
        std::vector<std::array<double, 6>> tracePoints; // x1, y1, z1, x2, y2, z2
        std::vector<bool> tips; // vero se non passante, falso se passante
        std::vector<double> lengths; // lunghezza della traccia
    };
}

