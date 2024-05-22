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
    DFN prova;
    prova.numFratture=8;
    prova.idFratture={0,1,2,3,4,5,6,7};
    prova.numVertici={4,3,4,3,3,4,3,5};
    prova.vertici={{Vector3d(2.5,2.,0.),Vector3d(-3.5,0.5,0.),Vector3d(-5.,4.5,0.),Vector3d(-1.,4.5,0.)},
                     {Vector3d(6.,1.,0.),Vector3d(6.,3.,1.),Vector3d(3.,1.5,7.)},
                     {Vector3d(5.,0.,0.),Vector3d(7.,0.,0.),Vector3d(7.,2.,0.),Vector3d(5.,2.,0.)},
                     {Vector3d(-4.,3.,0.2),Vector3d(-1.5,4.,2.),Vector3d(0.,0.,-1.)},
                     {Vector3d(0.,0.,-1),Vector3d(-4.,-3.,0.),Vector3d(-1.,-3.5,0.)},
                     {Vector3d(-4.,-3.,0.),Vector3d(-1.,-3.5,0.),Vector3d(-1.,-3.5,4.),Vector3d(-4.,-3.,3.)},
                     {Vector3d(1.5,-8.5,-3.),Vector3d(3.,-7.5,4.),Vector3d(7.,-3.5,0.)},
                     {Vector3d(4.5,-2.5,0.),Vector3d(3.5,-5.,0.),Vector3d(5.,-7.,0.),Vector3d(9.,-4.,0.),Vector3d(7.,-1.,0.)},
                     };
    memorizza_tracce(prova,tol);
    return 0;

}
    // Creazione di un'istanza di DFN con dati di esempio

// testare il primo output: OK
//     // Creazione di un'istanza della struttura DFN
//     LibraryDFN::DFN dfn;

//     // Popolamento dei dati della struttura (esempio di dati casuali)
//     dfn.numTracce = 2;
//     dfn.idTracce = {1, 2};
//     dfn.tracce = {{1, 2}, {2, 3}};
//     dfn.estremiTracce = {{0.0, 0.0, 0.0, 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0, 2.0, 2.0, 2.0}};

//     // Chiamata alla funzione per stampare le tracce su un file
//     printTraces(dfn, "tracce.txt");

//     std::cout << "Tracce stampate con successo." << std::endl;

//     // return 0;


    //testare secondo output: OK
//Creazione di un'istanza di DFN con dati di esempio
//     LibraryDFN::DFN dfn;
//     dfn.numFratture = 2;
//     dfn.idFratture = {1, 2};
//     dfn.numTracce = 4;
//     dfn.tracce = {{1, 2}, {2, 3}, {1, 3}, {2, 4}};
//     dfn.tips = {{{true, false}}, {{false, true}}, {{true, false}}, {{false, true}}};
//     dfn.lunghezze = {1.0, 2.0, 1.5, 1.8};

//     // Chiamata alla funzione per stampare le tracce per frattura
//     printTracesByFracture(dfn, "Output2.txt");
//     std::cout << "Seconda stampa terminata con successo." << std::endl;

//     return 0;
// }