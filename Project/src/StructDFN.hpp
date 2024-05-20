#pragma once

#include <vector>
#include "Eigen/Eigen"

namespace LibraryDFN{
struct DFN {

    //fratture
    unsigned int numFratture; //è il numero delle fratture, serve ad un bel niente (serve)
    std::vector<unsigned int> idFratture; //vettore degli id dei vettori (bisogna leggerli, non serve aggiungerne di nuovi dopo la fine della lettura da file)
    std::vector<unsigned int> numVertici; //vettore di numeri di vertici di ogni frattura (ordinati nel modo numVertici[i]= numero vertici della frattura idFratture[i]) (se non mi sbaglio serve a qualcosa, in lettura, non in scrittura)
    std::vector<std::vector<Eigen::Vector3d>> vertici; //vettore di matrici di dimensione 3x(numero vertici), ogni array individua le coordinate x, y, z del punto (ordinati nel modo vertici[i]= insieme delle coordinate dei vertici della frattura idFratture[i] in senso antiorario) (serve in lettura, non in scrittura)

    //tracce
    unsigned int numTracce; //è numero delle tracce (va stampato)
    std::vector<unsigned int> idTracce; //vettore degli id delle tracce (vanno letti)
    std::vector<std::array<unsigned int,2>> tracce; // idFratture1, idFratture2
    std::vector<std::array<double, 6>> estremiTracce; // x1, y1, z1, x2, y2, z2
    std::vector<std::array<bool,2>> tips; // vero se non passante, falso se passante
    std::vector<double> lunghezze; // lunghezza della traccia

};
}
