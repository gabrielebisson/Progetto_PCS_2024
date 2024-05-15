#include "Utils.hpp"
#include <iostream>

int main() {
    std::string filename = "DFN/FR3_data.txt";
    LibraryDFN::DFN dfn = LibraryDFN::readDFNFromFile(filename);

    if (dfn.numFratture == 0) {
        std::cerr << "Error: unable to read DFN from file" << std::endl;
        return 1;
    }

    std::cout << "Number of fractures: " << dfn.numFratture << std::endl;
    for (unsigned int i = 0; i < dfn.numFratture; ++i) {
        std::cout << "Fracture " << dfn.idFratture[i] << " with " << dfn.numVertici[i] << " vertices:" << std::endl;
        for (const auto& vertex : dfn.vertici[i]) {
            std::cout << "  (" << vertex.x() << ", " << vertex.y() << ", " << vertex.z() << ")" << std::endl;
        }
    }

    return 0;
}
