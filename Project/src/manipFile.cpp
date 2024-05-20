// Funzioni per manipolare i file di input e di output
#pragma once

#include "Utils.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

namespace LibraryDFN {

// ------------------- FUNZIONE LETTURA INPUT DA FILE -------------------------------

void readDFNFromFile(const std::string& filename, DFN& dfn)        //io farei che prendere il DFN in input e non lo creerei nella funzione, per risparmiare memoria, sennò poi gli tocca copiare tutto quando chiamate la funzione
{
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Errore: impossibile aprire il file " << filename << " ." << std::endl;
        return;
    }

    std::string line;

    // Leggere le fratture
    std::getline(file, line);
    std::getline(file, line);
    dfn.numFratture = std::stoi(line);


    // Riservare memoria per i vettori
    dfn.idFratture.reserve(dfn.numFratture);
    dfn.numVertici.reserve(dfn.numFratture);
    dfn.vertici.reserve(dfn.numFratture);

    // Leggere le fratture
    for (unsigned int i = 0; i < dfn.numFratture; ++i) {
        std::getline(file, line);
        std::getline(file, line);
        std::istringstream iss(line);
        unsigned int id;
        char sep;
        unsigned int numVertices;
        iss >> id >> sep >> numVertices;

        dfn.idFratture.push_back(id);                                                      //non avete riservato memoria per il vettore idfratture, così al programma tocca copiare ogni volta il vettore da qualche parte prima di poter aggiungere il nuovo id, sprecando memoria e tempo: fate un reserve all'inizio
        dfn.numVertici.push_back(numVertices);                                           //idem

        // Leggere i vertici
        std::vector<Eigen::Vector3d> vertices;
        vertices.reserve(numVertices);
        std::getline(file,line);
        std::istringstream iss2(line);

        for (int j = 0; j < 3; ++j) {
            if (!std::getline(file, line)){
                std::cerr << "Errore: impossibile leggere i vertici" << " ." << std::endl;
                break;
            }
            std::istringstream iss2(line);
            for (unsigned int k = 0; k < numVertices; ++k) {
                if (vertices.size() < numVertices)
                {
                    vertices.emplace_back(Eigen::Vector3d(0., 0., 0.));
                }
                iss2 >> vertices[k](j) >> sep;
                    }
                }

                dfn.vertici.push_back(vertices);                                            //di nuovo non avete riservato memoria per il vettore vertici...
    }

    file.close();
    return;
}

//--------------------- FUNZIONE STAMPA DELLE TRACCE ----------------------

void printTraces(const LibraryDFN::DFN& dfn, const std::string& filename)
{
    // Crea e apre un file di output
    std::ofstream outFile(filename);

    // Verifica se il file è stato aperto correttamente
    if (!outFile.is_open()) {
        std::cerr << "Errore: impossibile aprire il file " << filename << " ." << std::endl;
        return;
    }

    // Scrive il numero di tracce nel file
    outFile << "# Number of Traces" << dfn.numTracce << std::endl;

    // Scrive l'intestazione delle colonne
    outFile << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2"<< std::endl;

    // Stampa i valori voluti
    for (unsigned int i = 0; i < dfn.numTracce; ++i) {
        outFile << dfn.idTracce[i] << "; "
                << dfn.tracce[i][0] << "; " << dfn.tracce[i][1] << "; "
                << dfn.estremiTracce[i][0] << "; " << dfn.estremiTracce[i][1] << "; "
                << dfn.estremiTracce[i][2] << "; " << dfn.estremiTracce[i][3] << "; "
                << dfn.estremiTracce[i][4] << "; " << dfn.estremiTracce[i][5] << std::endl;
    }

    // Chiude il file
    outFile.close();
}


//----------------- FUNZIONE STAMPA TRACCE PASSANTI - NON PASSANTI -------------------------

void printTracesByFracture(const DFN& dfn, const std::string& filename)
{
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Errore nell'apertura del file: " << filename << " ." << std::endl;
        return;
    }

}
}
