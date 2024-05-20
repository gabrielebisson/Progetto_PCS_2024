#include "Utils.hpp"
#include "readFileTest.hpp"
#include <iostream>

//--------------------- TEST PER manipFile-------------------------
int main(){
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

    LibraryDFN::DFN dfn;
    dfn.numFratture = 2;
    dfn.idFratture = {1, 2};
    dfn.numVertici = {4, 3};
    dfn.vertici = {
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}},
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}}
    };
    dfn.numTracce = 4;
    dfn.idTracce = {1, 2, 3, 4};
    dfn.tracce = {{1, 2}, {2, 3}, {1, 3}, {2, 4}};
    dfn.estremiTracce = {{0.0, 0.0, 0.0, 1.0, 0.0, 0.0}, {1.0, 0.0, 0.0, 1.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0, 1.0, 0.0}, {1.0, 0.0, 0.0, 0.0, 1.0, 0.0}};
    dfn.tips = {true, false, true, false};
    dfn.lunghezze = {1.0, 1.0, 1.4142, 1.0}; // lunghezza della traccia

    // Chiamata alla funzione per stampare le tracce per frattura
    LibraryDFN::printTracesByFracture(dfn, "result.txt");

    return 0;
}

// ----------------------- G TEST PER readFile-------------------------
// int main(int argc, char **argv)
// {
//     testing::InitGoogleTest(&argc, argv);

//     return RUN_ALL_TESTS();
// }


