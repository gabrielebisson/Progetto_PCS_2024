#include "UCDUtilities.hpp"
#include "Utils.hpp"
#include "readFileTest.hpp"
#include <iostream>
#include <gtest/gtest.h>
#include <random>
#include "readFileTest.hpp"
#include "StructDFN.hpp"
#include "Eigen/Eigen"
#include "ExportDFN.cpp"


int main()
{

    //----------------------------- parte 1 ---------------------------------------

    double tol=__DBL_EPSILON__;
    LibraryDFN::DFN dfn;
    std::string nome_file="DFN/FR200_data.txt";
    readDFNFromFile(nome_file,dfn);
    memorizza_tracce(dfn,tol);

    Eigen::MatrixXd points;  // Matrice dei punti dei segmenti
    std::vector<Gedim::UCDProperty<double>> points_properties;  // Proprietà dei punti, se necessario
    std::vector<Gedim::UCDProperty<double>> segments_properties;  // Proprietà dei segmenti, se necessario
    Eigen::VectorXi materials;  // Materiali dei punti o dei segmenti, se necessario

    // cambio struttura dati dei vertici
    matriceVertici(dfn, points);

    // cambio struttura dati dei segmenti
    Eigen::MatrixXi segments;
    Eigen::MatrixXd punti;
    trasformaEstremiTracce(dfn.estremiTracce, punti); // metto solo i punti delle tracce
    creaMatriceSegments(dfn.estremiTracce, dfn.numTracce, punti, segments);

    // cambio struttura dati delle fratture
    std::vector<std::vector<unsigned int>> polygons_vertici;
    creaPolygonsVertices(dfn.vertici, points, polygons_vertici);

    // Imposta i materiali
    materials.setZero(); // Imposta tutti i materiali a 0

    Gedim::UCDUtilities ucdUtils;

    // esportazione dei vertici
    ucdUtils.ExportPoints("visualPara/FR3_VERTICI.ucd", points, points_properties, materials);
    //esportazione dei segmenti
    ucdUtils.ExportSegments("visualPara/FR3_TRACCE.ucd", punti, segments, points_properties, segments_properties, materials);
    //esportazione delle fratture
    ucdUtils.ExportPolygons("visualPara/FR3_FRATTURE.ucd", points, polygons_vertici);


//------------------------ parte 2 ------------------------------------

    sortTracesAndPrintByFracture(dfn,"TerrorTraxxPerFract.txt");
    definisci_mesh(dfn,tol);

    unsigned int maxTracceFratturaID = trovaFratturaMaxTracce(dfn);
    // std::cout << "La frattura con il maggior numero di tracce ha ID: " << maxTracceFratturaID << std::endl;

    // definizione mesh
    const PolygonalMesh& mesh = dfn.meshPoligonali[maxTracceFratturaID];

    // celle 0D
    Eigen::MatrixXd points0D = exportCell0DMesh(mesh);
    //std::cout << "Matrix:" << points0D << std::endl;
    Eigen::VectorXi materials0D(mesh.Cell0DMarkers.size());
    for (size_t i = 0; i < mesh.Cell0DMarkers.size(); ++i) {
        materials0D(i) = mesh.Cell0DMarkers[i];
    }


    // celle 1D
    Eigen::MatrixXi segments1D = convertToEigenMatrix(mesh.Cell1DVertices);
    Eigen::VectorXi materials1D(mesh.Cell1DMarkers.size());
    for (size_t i = 0; i < mesh.Cell1DMarkers.size(); ++i) {
        materials1D(i) = mesh.Cell1DMarkers[i];
    }

    Eigen::VectorXi materials2D(mesh.NumberCell2D);
    // Riempimento del vettore con valori incrementali
    for (unsigned int i = 0; i < mesh.NumberCell2D; ++i) {
        materials2D[i] = i + 1;
    }

    // Esporta i vertici celle 0D
    ucdUtils.ExportPoints("visualPara/FR200_cell0D.ucd", points0D, {}, materials0D);
    // Esporta i segmenti (Celle 1D)
    ucdUtils.ExportSegments("visualPara/FR200_cell1D.ucd", points0D, segments1D, {}, {}, materials1D);
    // esporta i poligoni (celle 2D)
    ucdUtils.ExportPolygons("visualPara/FR200_cell2D.ucd", points0D, mesh.Cell2DVertices, {}, {}, materials2D);

    return 0;

}
