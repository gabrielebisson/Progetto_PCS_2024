#include "Utils.hpp"
#include "readFileTest.hpp"
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
    std::string nome_file="DFN/FR200_data.txt";
    readDFNFromFile(nome_file,prova);
    memorizza_tracce(prova,tol);
    printTraces(prova,"risultato_tracce.txt");
    printTracesByFracture(prova,"risultato_tracce_ordinate.txt");
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
