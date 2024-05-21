#include "Utils.hpp"
#include "readFileTest.hpp"
#include <iostream>

//--------------------- TEST PER manipFile-------------------------
int main(){

// testiamo la lettura dei file: OK
    // std::string filename = "DFN/FR50_data.txt";
    // LibraryDFN::DFN dfn;
    // LibraryDFN::readDFNFromFile(filename, dfn);

    // if (dfn.numFratture == 0) {
    //     std::cerr << "Error: unable to read DFN from file" << std::endl;
    //     return 1;
    // }

    // std::cout << "Number of fractures: " << dfn.numFratture << std::endl;
    // for (unsigned int i = 0; i < dfn.numFratture; ++i) {
    //     std::cout << "Fracture " << dfn.idFratture[i] << " with " << dfn.numVertici[i] << " vertices:" << std::endl;
    //     for (const auto& vertex : dfn.vertici[i]) {
    //         std::cout << "  (" << vertex.x() << ", " << vertex.y() << ", " << vertex.z() << ")" << std::endl;
    //     }
    // }
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
    LibraryDFN::DFN dfn;
    dfn.numFratture = 2;
    dfn.idFratture = {1, 2};
    dfn.numTracce = 4;
    dfn.tracce = {{1, 2}, {2, 3}, {1, 3}, {2, 4}};
    dfn.tips = {{{true, false}}, {{false, true}}, {{true, false}}, {{false, true}}};
    dfn.lunghezze = {1.0, 2.0, 1.5, 1.8};

    // Chiamata alla funzione per stampare le tracce per frattura
    printTracesByFracture(dfn, "Output2.txt");
    std::cout << "Seconda stampa terminata con successo." << std::endl;

    return 0;
}

// ----------------------- G TEST PER readFile-------------------------
// int main(int argc, char **argv)
// {
//     testing::InitGoogleTest(&argc, argv);

//     return RUN_ALL_TESTS();
// }


