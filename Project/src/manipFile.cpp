// Funzioni per manipolare i file di input e di output

#include "Utils.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

namespace LibraryDFN {
DFN readDFNFromFile(const std::string& filename) {
    DFN dfn;
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: unable to open file " << filename << std::endl;
        return dfn;
    }

    std::string line;

    // Leggere il numero di fratture
    std::getline(file, line);
    std::getline(file, line);
    dfn.numFratture = std::stoi(line);
    //std::cout << dfn.numFratture<<std::endl;

    // Leggere le fratture
    for (unsigned int i = 0; i < dfn.numFratture; ++i) {
        std::getline(file, line);
        std::getline(file, line);
        std::istringstream iss(line);
        unsigned int id;
        char sep;
        unsigned int numVertices;
        iss >> id>> sep >> numVertices;
        // std::cout << id << "    " << numVertices << std::endl;
        dfn.idFratture.push_back(id);
        dfn.numVertici.push_back(numVertices);

        // Leggere i vertici
        std::vector<Eigen::Vector3d> vertices;
        vertices.reserve(numVertices);
            std::getline(file, line);
            std::istringstream iss2(line);
            for (int j = 0; j < 3; ++j) {
                if (!std::getline(file, line)) {
                    std::cerr << "Error: unable to read vertices" << std::endl;
                    break;
                }
                std::istringstream iss2(line);
                for (unsigned int k = 0; k < numVertices; ++k) {
                    if (vertices.size() < numVertices) {
                        vertices.emplace_back(Eigen::Vector3d(0., 0., 0.));
                    }
                    iss2 >> vertices[k](j) >> sep;
                }
            }

            dfn.vertici.push_back(vertices);
    }

    file.close();
    return dfn;
}
}
