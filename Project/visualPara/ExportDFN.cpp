#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "StructDFN.hpp"
#include "UCDUtilities.hpp"
#include <vector>
#include <Eigen/Dense>


//-------------------------------- TRASFORMA PARTE 1 ------------------------------------------

inline void matriceVertici(const LibraryDFN::DFN& dfn, Eigen::MatrixXd& points) {
    // Conta il numero totale di punti
    size_t totalPoints = 0;
    for(const auto& vec : dfn.vertici) {
        totalPoints += vec.size();
    }

    // Inizializza la matrice con il numero totale di punti
    points.resize(3, totalPoints);

    // Riempie la matrice con i dati da vertici
    size_t column = 0;
    for(const auto& vec : dfn.vertici) {
        for(const auto& point : vec) {
            points.col(column++) = point;
        }
    }

}

inline void trasformaEstremiTracce(const std::vector<std::array<Eigen::Vector3d, 2>>& estremiTracce, Eigen::MatrixXd& points) {
    // Calcola il numero totale di punti (righe della matrice)
    size_t totalPoints = estremiTracce.size() * 2; // ogni array ha due elementi

    // Inizializza la matrice con il numero totale di punti (3 righe per ogni punto)
    points.resize(3, totalPoints);

    // Riempie la matrice con i dati da estremiTracce
    size_t column = 0;
    for(const auto& arr : estremiTracce) {
        points.col(column++) = arr[0]; // primo elemento dell'array
        points.col(column++) = arr[1]; // secondo elemento dell'array
    }

}


inline void creaMatriceSegments(const std::vector<std::array<Eigen::Vector3d, 2>>& estremiTracce, unsigned int numTracce,
                                const Eigen::MatrixXd& points,
                                Eigen::MatrixXi& segments)
{
    // Inizializza la matrice segments con dimensioni 2 x numTracce
    segments.resize(2, numTracce);

    // Riempie la matrice segments con gli indici dei vertici iniziali e finali
    for (unsigned int i = 0; i < numTracce; ++i) {
        // Trova gli indici dei vertici iniziale e finale nella matrice points
        Eigen::Vector3d verticeIniziale = estremiTracce[i][0];
        Eigen::Vector3d verticeFinale = estremiTracce[i][1];

        // Cerca gli indici dei vertici iniziale e finale nella matrice points
        int indiceIniziale = -1;
        int indiceFinale = -1;

        for (int j = 0; j < points.cols(); ++j) {
            if (points.col(j) == verticeIniziale) {
                indiceIniziale = j;
            }
            if (points.col(j) == verticeFinale) {
                indiceFinale = j;
            }
            if (indiceIniziale != -1 && indiceFinale != -1) {
                break;  // Esci dal ciclo se hai trovato entrambi gli indici
            }
        }

        // Assegna gli indici alla matrice segments
        segments(0, i) = indiceIniziale;
        segments(1, i) = indiceFinale;
    }

}


inline void creaPolygonsVertices(const std::vector<std::vector<Eigen::Vector3d>>& vertici,
                                 const Eigen::MatrixXd& points,
                                 std::vector<std::vector<unsigned int>>& polygons_vertices) {
    // Numero delle fratture
    unsigned int numFratture = vertici.size();

    // Inizializza il vettore polygons_vertices con la dimensione di numFratture
    polygons_vertices.resize(numFratture);

    // Popola polygons_vertices con gli ID dei vertici
    for (unsigned int i = 0; i < numFratture; ++i) {
        for (const auto& vertex : vertici[i]) {
            // Trova l'ID del vertice confrontando con i punti nella matrice points
            for (unsigned int j = 0; j < points.cols(); ++j) {
                if (vertex == points.col(j)) {
                    polygons_vertices[i].push_back(j);
                    break;
                }
            }
        }
    }
}


//-------------------------------- TRASFORMA PARTE 2 ------------------------------------------

// Trova la frattura con il maggior numero di tracce in modo efficiente
inline unsigned int trovaFratturaMaxTracce(const LibraryDFN::DFN& dfn) {
    unsigned int maxTraces = 0;
    unsigned int maxTracesFractureId = 0;

    for (unsigned int i = 0; i < dfn.numFratture; ++i) {
        unsigned int numTracesForFracture = 0;
        for (unsigned int j = 0; j < dfn.numTracce; ++j) {
            if (dfn.tracce[j][0] == dfn.idFratture[i] || dfn.tracce[j][1] == dfn.idFratture[i]) {
                numTracesForFracture++;
            }
        }
        if (numTracesForFracture > maxTraces) {
            maxTraces = numTracesForFracture;
            maxTracesFractureId = dfn.idFratture[i];
        }

    }

    //std::cout << "Max traces: " << maxTraces << std::endl;

    return maxTracesFractureId;
}


// Ottieni la mesh della frattura data il suo ID
inline const PolygonalLibrary::PolygonalMesh& getFractureMesh(const LibraryDFN::DFN& dfn, unsigned int fractureID) {
    for (unsigned int i = 0; i < dfn.idFratture.size(); ++i) {
        if (dfn.idFratture[i] == fractureID) {
            return dfn.meshPoligonali[i];
        }
    }
    throw std::runtime_error("Fracture ID not found");
}


inline Eigen::MatrixXd exportCell0DMesh(const PolygonalLibrary::PolygonalMesh& mesh ) {

    Eigen::MatrixXd points0D;

    // Determinare il numero di punti
    int num_points = mesh.Cell0DCoordinates.size();

    // Creare una matrice Eigen::MatrixXd con dimensioni appropriate
    points0D.resize(3, num_points);

    // Riempire la matrice con i dati dal vettore
    for (int i = 0; i < num_points; ++i) {
        points0D.col(i) = mesh.Cell0DCoordinates[i];
    }

    // Stampa la matrice per verificarne il contenuto
    // std::cout << "Matrix:" << points0D << std::endl;

    return points0D;

}


inline Eigen::MatrixXi convertToEigenMatrix(const std::vector<std::array<unsigned int,2>>& Cell1DVertices) {
    std::size_t numVertices = Cell1DVertices.size();

    // Creiamo una matrice con 2 righe e numVertices colonne
    Eigen::MatrixXi lines(2, numVertices);

    for (std::size_t i = 0; i < numVertices; ++i) {
        lines(0, i) = static_cast<int>(Cell1DVertices[i][0]); // ID_inizio
        lines(1, i) = static_cast<int>(Cell1DVertices[i][1]); // ID_fine
    }

    return lines;
}








