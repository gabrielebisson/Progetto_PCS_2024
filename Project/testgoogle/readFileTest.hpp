#pragma once

#include <gtest/gtest.h>
#include <fstream>
#include "StructDFN.hpp"
#include "PolygonalMesh.hpp"
#include "Utils.hpp"
#include <Eigen/Eigen>
#include <vector>
#include <array>
#include <tuple>

// Test per la lettura del DFN da un file
TEST(ReadDFNFromFileTest, FileInput) {
    LibraryDFN::DFN dfn;

    // Chiamare la funzione per leggere il DFN dal file
    LibraryDFN::readDFNFromFile("DFN/FR3_data.txt", dfn);
    // Verificare il numero di fratture
    EXPECT_EQ(dfn.numFratture, 3);

    // Verificare l'ID delle fratture e il numero di vertici
    ASSERT_EQ(dfn.idFratture.size(), 3);
    ASSERT_EQ(dfn.numVertici.size(), 3);
    ASSERT_EQ(dfn.vertici.size(), 3);

    // Prima frattura
    EXPECT_EQ(dfn.idFratture[0], 0);
    EXPECT_EQ(dfn.numVertici[0], 4);
    ASSERT_EQ(dfn.vertici[0].size(), 4);
    EXPECT_EQ(dfn.vertici[0][0], Eigen::Vector3d(0.0, 0.0, 0.0));
    EXPECT_EQ(dfn.vertici[0][1], Eigen::Vector3d(1.0, 0.0, 0.0));
    EXPECT_EQ(dfn.vertici[0][2], Eigen::Vector3d(1.0, 1.0, 0.0));
    EXPECT_EQ(dfn.vertici[0][3], Eigen::Vector3d(0.0, 1.0, 0.0));

    // Seconda frattura
    EXPECT_EQ(dfn.idFratture[1], 1);
    EXPECT_EQ(dfn.numVertici[1], 4);
    ASSERT_EQ(dfn.vertici[1].size(), 4);
    EXPECT_EQ(dfn.vertici[1][0], Eigen::Vector3d(0.8, 0.0, -0.1));
    EXPECT_EQ(dfn.vertici[1][1], Eigen::Vector3d(0.8, 0.0, 0.29999999999999999));
    EXPECT_EQ(dfn.vertici[1][2], Eigen::Vector3d(0.8, 1.0, 0.29999999999999999));
    EXPECT_EQ(dfn.vertici[1][3], Eigen::Vector3d(0.8, 1.0, -0.1));

    // Terza frattura
    EXPECT_EQ(dfn.idFratture[2], 2);
    EXPECT_EQ(dfn.numVertici[2], 4);
    ASSERT_EQ(dfn.vertici[2].size(), 4);
    EXPECT_EQ(dfn.vertici[2][0], Eigen::Vector3d(-0.23777799999999999, 0.5, -0.34444000000000002));
    EXPECT_EQ(dfn.vertici[2][1], Eigen::Vector3d(0.31618370000000001, 0.5, -0.34444000000000002));
    EXPECT_EQ(dfn.vertici[2][2], Eigen::Vector3d(0.31618370000000001, 0.5, 0.45283889999999999));
    EXPECT_EQ(dfn.vertici[2][3], Eigen::Vector3d(-0.23777799999999999, 0.5, 0.45283889999999999));
}



// Test per la funzione printTraces
TEST(PrintTracesTest, FileOutput) {
    LibraryDFN::DFN dfn;

    // Dati di test
    dfn.numTracce = 2;
    dfn.idTracce = {0, 1};
    dfn.tracce = {{0, 1}, {1, 2}};
    dfn.estremiTracce = {
        {{Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(1.0, 1.0, 1.0)}},
        {{Eigen::Vector3d(1.0, 1.0, 1.0), Eigen::Vector3d(2.0, 2.0, 2.0)}}
    };

    std::string tempFile = "temp_output.txt"; // file temporaneo
    printTraces(dfn, tempFile);

    std::ifstream inFile(tempFile); // leggo file
    ASSERT_TRUE(inFile.is_open());

    std::string line;
    std::vector<std::string> lines;
    while (std::getline(inFile, line)) {
        lines.push_back(line);
    }
    inFile.close();

    std::remove(tempFile.c_str());

    // Verifico output atteso
    ASSERT_EQ(lines.size(), 5);
    EXPECT_EQ(lines[0], "# Number of Traces");
    EXPECT_EQ(lines[1], "2");
    EXPECT_EQ(lines[2], "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2");
    EXPECT_EQ(lines[3], "0; 0; 1; 0; 0; 0; 1; 1; 1");
    EXPECT_EQ(lines[4], "1; 1; 2; 1; 1; 1; 2; 2; 2");
}

// Test per la funzione printTracesByFracture
TEST(PrintTracesByFractureTest, Output) {
    LibraryDFN::DFN dfn;

    // Dati di test
    dfn.numFratture = 3;
    dfn.idFratture = {0, 1, 2};
    dfn.numTracce = 3;
    dfn.idTracce = {0, 1, 2};
    dfn.tracce = {{0, 1}, {1, 0}, {0, 1}};
    dfn.tips = {{false, true}, {true, false}, {false, false}};
    dfn.lunghezze = {10.0, 20.0, 15.0};

    std::string tempFile = "temp_fracture_output.txt";
    printTracesByFracture(dfn, tempFile);

    std::ifstream inFile(tempFile);
    ASSERT_TRUE(inFile.is_open());

    std::string line;
    std::vector<std::string> lines;
    while (std::getline(inFile, line)) {
        lines.push_back(line);
    }
    inFile.close();

    std::remove(tempFile.c_str());

    ASSERT_EQ(lines.size(), 15);

    // prima frattura
    EXPECT_EQ(lines[0], "# FractureId; NumTraces");
    EXPECT_EQ(lines[1], "0; 3");
    EXPECT_EQ(lines[2], "# TraceId; Tips; Length");
    EXPECT_EQ(lines[3], "1; 1; 20");
    EXPECT_EQ(lines[4], "2; 0; 15");
    EXPECT_EQ(lines[5], "0; 0; 10");

    // seconda frattura
    EXPECT_EQ(lines[6], "# FractureId; NumTraces");
    EXPECT_EQ(lines[7], "1; 3");
    EXPECT_EQ(lines[8], "# TraceId; Tips; Length");
    EXPECT_EQ(lines[9], "0; 1; 10");
    EXPECT_EQ(lines[10], "1; 0; 20");
    EXPECT_EQ(lines[11], "2; 0; 15");

    // terza frattura
    EXPECT_EQ(lines[12], "# FractureId; NumTraces");
    EXPECT_EQ(lines[13], "2; 0");
    EXPECT_EQ(lines[14], "# TraceId; Tips; Length");
}



// Test per funzione che triangola la frattura
TEST(TriangolaFratturaTest, Frattura) {
    LibraryDFN::DFN dfn;
    dfn.numVertici = {5};  // 1 frattura con 5 vertici

    std::vector<std::array<unsigned int, 3>> expected = {{0, 1, 2}, {0, 2, 3}, {0, 3, 4}};
    auto result = triangola_frattura(dfn, 0);
    EXPECT_EQ(result, expected);
}



// Test per funzione che trova il versore perpendicolare al piano che contiene il poligono
TEST(VersoreNormaleTest, Poligono) {
    std::vector<Eigen::Vector3d> poligono = {
        Eigen::Vector3d(0, 0, 0),
        Eigen::Vector3d(1, 0, 0),
        Eigen::Vector3d(1, 1, 0),
        Eigen::Vector3d(0, 1, 0)
    };

    Eigen::Vector3d expected_normal(0, 0, 1);
    Eigen::Vector3d normal = LibraryDFN::versore_normale(poligono);

    EXPECT_NEAR(normal[0], expected_normal[0], 1e-6);
    EXPECT_NEAR(normal[1], expected_normal[1], 1e-6);
    EXPECT_NEAR(normal[2], expected_normal[2], 1e-6);
}



// Test funzione scarta fratture
// Test per il caso in cui le fratture non si intersecano
TEST(ScartaFrattureTest, NoIntersection) {
    LibraryDFN::DFN dfn;
    dfn.numFratture = 2;
    dfn.vertici = {
        {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(1, 1, 0), Eigen::Vector3d(0, 1, 0)},
        {Eigen::Vector3d(2, 2, 0), Eigen::Vector3d(3, 2, 0), Eigen::Vector3d(3, 3, 0), Eigen::Vector3d(2, 3, 0)}
    };
    dfn.numVertici = {4, 4};
    dfn.versori = {Eigen::Vector3d(0, 0, 1), Eigen::Vector3d(0, 0, 1)};
    double tol = 0.1;

    auto result = LibraryDFN::scarta_fratture(dfn, tol);
    EXPECT_TRUE(result.empty());
}

//Test per il caso in cui le fratture si intersecano
TEST(ScartaFrattureTest, Intersection) {
    LibraryDFN::DFN dfn;
    dfn.numFratture = 2;
    dfn.vertici = {
        {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(2, 0, 0), Eigen::Vector3d(2, 2, 0), Eigen::Vector3d(0, 2, 0)},
        {Eigen::Vector3d(1, 1, -1), Eigen::Vector3d(3, 1, -1), Eigen::Vector3d(3, 1, 1), Eigen::Vector3d(1, 1, 1)}
    };
    dfn.numVertici = {4, 4};
    dfn.versori = {Eigen::Vector3d(0, 0, 1), Eigen::Vector3d(0, 1, 0)};
    double tol = 0.1;

    auto result = LibraryDFN::scarta_fratture(dfn, tol);

    ASSERT_EQ(result.size(), 1);
    EXPECT_EQ(result[0][0], 0);
    EXPECT_EQ(result[0][1], 1);
}

// Test per il caso in cui le fratture si toccano ma non si intersecano
TEST(ScartaFrattureTest, Touching) {
    LibraryDFN::DFN dfn;
    dfn.numFratture = 2;
    dfn.vertici = {
        {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(1, 1, 0), Eigen::Vector3d(0, 1, 0)},
        {Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(2, 0, 0), Eigen::Vector3d(2, 1, 0), Eigen::Vector3d(1, 1, 0)}
    };
    dfn.numVertici = {4, 4};
    dfn.versori = {Eigen::Vector3d(0, 0, 1), Eigen::Vector3d(0, 0, 1)};
    double tol = 0.1;

    auto result = LibraryDFN::scarta_fratture(dfn, tol);
    EXPECT_TRUE(result.empty());
}



// Test per funzione che memorizza le tracce
TEST(MemorizzaTracceTest, DFN) {
    LibraryDFN::DFN dfn;

    // Definizione delle fratture e dei vertici
    dfn.numFratture = 3;
    dfn.idFratture = {0, 1, 2};
    dfn.numVertici = {4, 4, 4};
    dfn.vertici = {
        {
            Eigen::Vector3d(0.0, 0.0, 0.0),
            Eigen::Vector3d(1.0, 0.0, 0.0),
            Eigen::Vector3d(1.0, 1.0, 0.0),
            Eigen::Vector3d(0.0, 1.0, 0.0)
        },
        {
            Eigen::Vector3d(0.8, 0.0, -0.1),
            Eigen::Vector3d(0.8, 0.0, 0.29999999999999999),
            Eigen::Vector3d(0.8, 1.0, 0.29999999999999999),
            Eigen::Vector3d(0.8, 1.0, -0.1)
        },
        {
            Eigen::Vector3d(-0.23777799999999999, 0.5, -0.34444000000000002),
            Eigen::Vector3d(0.31618370000000001, 0.5, -0.34444000000000002),
            Eigen::Vector3d(0.31618370000000001, 0.5, 0.45283889999999999),
            Eigen::Vector3d(-0.23777799999999999, 0.5, 0.45283889999999999)
        }
    };

    memorizza_tracce(dfn, 1e-6);

    // Prima frattura
    EXPECT_EQ(dfn.idFratture[0], 0);
    EXPECT_EQ(dfn.numVertici[0], 4);
    ASSERT_EQ(dfn.vertici[0].size(), 4);
    EXPECT_EQ(dfn.vertici[0][0], Eigen::Vector3d(0.0, 0.0, 0.0));
    EXPECT_EQ(dfn.vertici[0][1], Eigen::Vector3d(1.0, 0.0, 0.0));
    EXPECT_EQ(dfn.vertici[0][2], Eigen::Vector3d(1.0, 1.0, 0.0));
    EXPECT_EQ(dfn.vertici[0][3], Eigen::Vector3d(0.0, 1.0, 0.0));

    // Seconda frattura
    EXPECT_EQ(dfn.idFratture[1], 1);
    EXPECT_EQ(dfn.numVertici[1], 4);
    ASSERT_EQ(dfn.vertici[1].size(), 4);
    EXPECT_EQ(dfn.vertici[1][0], Eigen::Vector3d(0.8, 0.0, -0.1));
    EXPECT_EQ(dfn.vertici[1][1], Eigen::Vector3d(0.8, 0.0, 0.29999999999999999));
    EXPECT_EQ(dfn.vertici[1][2], Eigen::Vector3d(0.8, 1.0, 0.29999999999999999));
    EXPECT_EQ(dfn.vertici[1][3], Eigen::Vector3d(0.8, 1.0, -0.1));

    // Terza frattura
    EXPECT_EQ(dfn.idFratture[2], 2);
    EXPECT_EQ(dfn.numVertici[2], 4);
    ASSERT_EQ(dfn.vertici[2].size(), 4);
    EXPECT_EQ(dfn.vertici[2][0], Eigen::Vector3d(-0.23777799999999999, 0.5, -0.34444000000000002));
    EXPECT_EQ(dfn.vertici[2][1], Eigen::Vector3d(0.31618370000000001, 0.5, -0.34444000000000002));
    EXPECT_EQ(dfn.vertici[2][2], Eigen::Vector3d(0.31618370000000001, 0.5, 0.45283889999999999));
    EXPECT_EQ(dfn.vertici[2][3], Eigen::Vector3d(-0.23777799999999999, 0.5, 0.45283889999999999));
}



// Test per la funzione che trova l'intersezione tra 2 segmenti
// Segmenti che si intersecano
TEST(IntersecaSegmentiTest, Intersection) {
    // Definizione dei punti dei segmenti A e B
    Eigen::Vector3d A1(0, 0, 0);
    Eigen::Vector3d A2(1, 1, 1);
    Eigen::Vector3d B1(0, 1, 0);
    Eigen::Vector3d B2(1, 0, 1);
    double tol = 1e-6;

    // Chiamata della funzione per verificare l'intersezione
    auto [point, beta] = LibraryDFN::interseca_segmenti(A1, A2, B1, B2);

    // Verifica che i segmenti si intersecano e che il punto di intersezione è corretto
    EXPECT_NEAR(point.x(), 0.5, tol);
    EXPECT_NEAR(point.y(), 0.5, tol);
    EXPECT_NEAR(point.z(), 0.5, tol);
    EXPECT_NEAR(beta, 0.5, tol); // Verifica se beta è all'interno dell'intervallo
}

// Segmenti che non si intersecano
TEST(IntersecaSegmentiTest, NoIntersection) {
    // Definizione dei punti dei segmenti A e B
    Eigen::Vector3d A1(0, 0, 0);
    Eigen::Vector3d A2(1, 1, 0);
    Eigen::Vector3d B1(2, 2, 0);
    Eigen::Vector3d B2(3, 3, 0);
    double tol = 1e-6;

    // Chiamata della funzione per verificare l'intersezione
    auto [point, beta] = LibraryDFN::interseca_segmenti(A1, A2, B1, B2);

    // Verifica che i segmenti non si intersecano
    EXPECT_TRUE(beta < 0 || beta > 1);
}



// Test per la funzione che aggiorna la mesh
TEST(AggiornaMeshTest, BaseCase) {
    // Imposta dati della mesh
    PolygonalMesh mesh;
    double tol = std::numeric_limits<double>::epsilon();
    mesh.NumberCell2D = 3;
    mesh.NumberCell1D = 5;
    mesh.NumberCell0D = 6;
    mesh.Cell2DId = {0, 1, 2};
    mesh.Cell1DId = {0, 1, 2, 3, 4};
    mesh.Cell0DId = {0, 1, 2, 3, 4, 5};
    mesh.Cell2DVertices = {{0, 1, 2, 3}, {2, 3, 4}, {3, 4, 5}};
    mesh.Cell2DEdges = {{0, 1, 2, 3}, {2, 3, 4}, {3, 4, 0}};
    mesh.Cell1DVertices = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}};
    mesh.Cell1DMarkers = {0, 1, 1, 0, 0};
    mesh.Cell0DCoordinates = {
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(1.0, 0.0, 0.0),
        Eigen::Vector3d(1.0, 1.0, 0.0),
        Eigen::Vector3d(0.0, 1.0, 0.0),
        Eigen::Vector3d(-1.0, 1.0, 0.0),
        Eigen::Vector3d(-1.0, 0.0, 0.0)
    };

    // Parametri di test per aggiorna_mesh
    std::vector<unsigned int> poligoni = {0};
    std::vector<std::array<unsigned int, 2>> lati_coinvolti = {{0, 1}, {1, 2}};
    Eigen::Vector3d est1(0.5, 0.0, 0.0);
    Eigen::Vector3d est2(0.5, 0.0, 0.0);

    // Imposta dati DFN
    LibraryDFN::DFN dfn;
    dfn.numFratture = 3;
    dfn.idFratture = {0, 1, 2};
    dfn.numVertici = {{7}, {7}, {3}};
    dfn.vertici = {
        {Eigen::Vector3d(5.0, 2.0, 0.0), Eigen::Vector3d(11.0, 0.0, 0.0), Eigen::Vector3d(20.0, 0.0, 0.0), Eigen::Vector3d(22.0, 8.0, 0.0), Eigen::Vector3d(14.0, 16.0, 0.0), Eigen::Vector3d(2.0, 13.0, 0.0), Eigen::Vector3d(0.0, 7.0, 0.0)},
        {Eigen::Vector3d(-19.0, 0.0, 0.0), Eigen::Vector3d(-7.0, 0.0, 0.0), Eigen::Vector3d(-5.0, 3.0, 0.0), Eigen::Vector3d(-7.0, 6.0, 0.0), Eigen::Vector3d(-10.0, 10.0, 0.0), Eigen::Vector3d(-13.0, 11.0, 0.0), Eigen::Vector3d(-19.0, 8.0, 0.0)},
        {Eigen::Vector3d(4.0, -8.0, 0.0), Eigen::Vector3d(12.0, -8.0, 0.0), Eigen::Vector3d(8.0, -1.0, 0.0)}
    };
    dfn.versori = {Eigen::Vector3d(0.0, 0.0, 1.0), Eigen::Vector3d(0.0, 0.0, 1.0), Eigen::Vector3d(0.0, 0.0, 1.0)};
    dfn.idTracce = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
    dfn.estremiTracce = {
        {Eigen::Vector3d(1.0, 10.0, 0.0), Eigen::Vector3d(17.0, 13.0, 0.0)},
        {Eigen::Vector3d(6.0, 14.0, 0.0), Eigen::Vector3d(8.0, 1.0, 0.0)},
        {Eigen::Vector3d(19.0, 11.0, 0.0), Eigen::Vector3d(19.0, 0.0, 0.0)},
        {Eigen::Vector3d(10.0, 13.0, 0.0), Eigen::Vector3d(6.0, 6.0, 0.0)},
        {Eigen::Vector3d(15.0, 9.0, 0.0), Eigen::Vector3d(19.0, 11.0, 0.0)},
        {Eigen::Vector3d(15.0, 6.0, 0.0), Eigen::Vector3d(19.0, 5.0, 0.0)},
        {Eigen::Vector3d(11.0, 7.0, 0.0), Eigen::Vector3d(12.0, 5.0, 0.0)},
        {Eigen::Vector3d(-13.0, 0.0, 0.0), Eigen::Vector3d(-13.0, 11.0, 0.0)},
        {Eigen::Vector3d(-18.0, 6.0, 0.0), Eigen::Vector3d(-10.0, 6.0, 0.0)},
        {Eigen::Vector3d(-11.0, 5.0, 0.0), Eigen::Vector3d(-10.0, 1.0, 0.0)},
        {Eigen::Vector3d(-13.0, 6.0, 0.0), Eigen::Vector3d(-15.0, 9.0, 0.0)},
        {Eigen::Vector3d(-17.0, 2.0, 0.0), Eigen::Vector3d(-15.0, 4.0, 0.0)},
        {Eigen::Vector3d(-11.5, 8.0, 0.0), Eigen::Vector3d(-13.0, 6.0, 0.0)},
        {Eigen::Vector3d(-7.0, 2.0, 0.0), Eigen::Vector3d(-7.0, 3.0, 0.0)},
        {Eigen::Vector3d(8.0, -5.0, 0.0), Eigen::Vector3d(12.0, -8.0, 0.0)},
        {Eigen::Vector3d(8.0, -3.5, 0.0), Eigen::Vector3d(8.0, -1.0, 0.0)},
        {Eigen::Vector3d(6.0, -6.5, 0.0), Eigen::Vector3d(4.0, -8.0, 0.0)},
        {Eigen::Vector3d(8.0, -6.0, 0.0), Eigen::Vector3d(8.0, -4.0, 0.0)},
        {Eigen::Vector3d(6.0, 14.0, 0.0), Eigen::Vector3d(10.0, 15.0, 0.0)}
    };
    dfn.traccePassanti = {{0, 1, 2}, {7}, {}};
    dfn.tracceNonPassanti = {{3, 4, 5, 6, 18}, {8, 9, 10, 11, 12, 13}, {14, 15, 16, 17}};

    // Chiamata alla funzione da testare
    LibraryDFN::definisci_mesh(dfn, tol);
    LibraryDFN::aggiorna_mesh(mesh, poligoni, lati_coinvolti, est1, est2, tol);

    // Verifiche sui risultati attesi
    ASSERT_EQ(mesh.NumberCell2D, 4);
    ASSERT_EQ(mesh.NumberCell1D, 7);
    ASSERT_EQ(mesh.NumberCell0D, 7);

    // Verifica la struttura delle celle 2D
    ASSERT_EQ(mesh.Cell2DVertices[0].size(), 5);
    ASSERT_EQ(mesh.Cell2DVertices[1].size(), 3);
    ASSERT_EQ(mesh.Cell2DVertices[2].size(), 4);
    ASSERT_EQ(mesh.Cell2DVertices[3].size(), 3);

    // Verifica la struttura delle celle 1D
    ASSERT_EQ(mesh.Cell1DVertices.size(), 7);
    ASSERT_EQ(mesh.Cell1DVertices[5], (std::array<unsigned int, 2>{6, 1}));
    ASSERT_EQ(mesh.Cell1DVertices[6], (std::array<unsigned int, 2>{6, 1}));

    // Verifica markers
    ASSERT_EQ(mesh.Cell1DMarkers.size(), 7);
    ASSERT_EQ(mesh.Cell1DMarkers[5], 0);
    ASSERT_EQ(mesh.Cell1DMarkers[6], 0);

}



// Test per la funzione DefinisciMesh
TEST(DefinisciMeshTest, BaseCase) {
    double tol = std::numeric_limits<double>::epsilon();

    // Definizione dei dati per DFN
    LibraryDFN::DFN dfn;
    dfn.numFratture = 3;
    dfn.idFratture = {0, 1, 2};
    dfn.numVertici = {{7}, {7}, {3}};
    dfn.vertici = {
        {
            Eigen::Vector3d(5., 2., 0.),
            Eigen::Vector3d(11., 0., 0.),
            Eigen::Vector3d(20., 0., 0.),
            Eigen::Vector3d(22., 8., 0.),
            Eigen::Vector3d(14., 16., 0.),
            Eigen::Vector3d(2., 13., 0.),
            Eigen::Vector3d(0., 7., 0.)
        },
        {
            Eigen::Vector3d(-19., 0, 0.),
            Eigen::Vector3d(-7., 0., 0.),
            Eigen::Vector3d(-5., 3., 0.),
            Eigen::Vector3d(-7., 6., 0.),
            Eigen::Vector3d(-10., 10., 0.),
            Eigen::Vector3d(-13., 11., 0.),
            Eigen::Vector3d(-19., 8., 0.)
        },
        {
            Eigen::Vector3d(4., -8., 0.),
            Eigen::Vector3d(12., -8., 0.),
            Eigen::Vector3d(8., -1., 0.)
        }
    };
    dfn.versori = {Eigen::Vector3d(0., 0., 1.), Eigen::Vector3d(0., 0., 1.), Eigen::Vector3d(0., 0., 1.)};
    dfn.idTracce = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
    dfn.estremiTracce = {
        {Eigen::Vector3d(1., 10., 0.), Eigen::Vector3d(17., 13., 0.)},
        {Eigen::Vector3d(6., 14., 0.), Eigen::Vector3d(8., 1., 0.)},
        {Eigen::Vector3d(19., 11., 0.), Eigen::Vector3d(19., 0., 0.)},
        {Eigen::Vector3d(10., 13., 0.), Eigen::Vector3d(6., 6., 0.)},
        {Eigen::Vector3d(15., 9., 0.), Eigen::Vector3d(19., 11., 0.)},
        {Eigen::Vector3d(15., 6., 0.), Eigen::Vector3d(19., 5., 0.)},
        {Eigen::Vector3d(11., 7., 0.), Eigen::Vector3d(12., 5., 0.)},
        {Eigen::Vector3d(-13., 0., 0.), Eigen::Vector3d(-13., 11., 0.)},
        {Eigen::Vector3d(-18., 6., 0.), Eigen::Vector3d(-10., 6., 0.)},
        {Eigen::Vector3d(-11., 5., 0.), Eigen::Vector3d(-10., 1., 0.)},
        {Eigen::Vector3d(-13., 6., 0.), Eigen::Vector3d(-15., 9., 0.)},
        {Eigen::Vector3d(-17., 2., 0.), Eigen::Vector3d(-15., 4., 0.)},
        {Eigen::Vector3d(-11.5, 8., 0.), Eigen::Vector3d(-13., 6., 0.)},
        {Eigen::Vector3d(-7., 2., 0.), Eigen::Vector3d(-7., 3., 0.)},
        {Eigen::Vector3d(8., -5., 0.), Eigen::Vector3d(12., -8., 0.)},
        {Eigen::Vector3d(8., -3.5, 0.), Eigen::Vector3d(8., -1., 0.)},
        {Eigen::Vector3d(6., -6.5, 0.), Eigen::Vector3d(4., -8., 0.)},
        {Eigen::Vector3d(8., -6., 0.), Eigen::Vector3d(8., -4., 0.)},
        {Eigen::Vector3d(6., 14., 0.), Eigen::Vector3d(10., 15., 0.)}
    };
    dfn.traccePassanti = {{0, 1, 2}, {7}, {}};
    dfn.tracceNonPassanti = {{3, 4, 5, 6, 18}, {8, 9, 10, 11, 12, 13}, {14, 15, 16, 17}};

    // Chiamata alla funzione da testare
    LibraryDFN::definisci_mesh(dfn, tol);

    // Verifiche sui risultati attesi
    // Numero di fratture
    ASSERT_EQ(dfn.numFratture, 3);

    // Numero di vertici per frattura
    ASSERT_EQ(dfn.numVertici[0], 7);
    ASSERT_EQ(dfn.numVertici[1], 7);
    ASSERT_EQ(dfn.numVertici[2], 3);

    for (size_t i = 0; i < dfn.vertici.size(); ++i) {
        for (size_t j = 0; j < dfn.vertici[i].size(); ++j) {
            ASSERT_TRUE(dfn.vertici[i][j].isApprox(expectedVertices[i][j], tol));
        }
    }

    // Verificare le tracce passanti e non passanti
    ASSERT_EQ(dfn.traccePassanti[0].size(), 3);
    ASSERT_EQ(dfn.traccePassanti[1].size(), 1);
    ASSERT_EQ(dfn.traccePassanti[2].size(), 0);

    ASSERT_EQ(dfn.tracceNonPassanti[0].size(), 5);
    ASSERT_EQ(dfn.tracceNonPassanti[1].size(), 6);
    ASSERT_EQ(dfn.tracceNonPassanti[2].size(), 4);
}


// Test per nuovo_poligono
// caso normale
TEST(NuovoPoligonoTest, CasoNormale) {
    // Inizializza PolygonalMesh
    PolygonalLibrary::PolygonalMesh mesh;
    mesh.NumberCell0D = 4;
    mesh.Cell0DId = {0, 1, 2, 3};
    mesh.Cell0DCoordinates = {Eigen::Vector3d(0,0,0), Eigen::Vector3d(1,0,0), Eigen::Vector3d(1,1,0), Eigen::Vector3d(0,1,0)};
    mesh.Cell0DMarkers = {0, 0, 0, 0};

    mesh.NumberCell1D = 4;
    mesh.Cell1DId = {0, 1, 2, 3};
    mesh.Cell1DVertices = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    mesh.Cell1DMarkers = {0, 0, 0, 0};

    mesh.NumberCell2D = 1;
    mesh.Cell2DId = {0};
    mesh.Cell2DVertices = {{0, 1, 2, 3}};
    mesh.Cell2DEdges = {{0, 1, 2, 3}};

    // Inizializza la mappa dei vecchi lati ai nuovi lati
    std::map<unsigned int, std::vector<unsigned int>> mappa_vecchi_lati_nuovi_lati = {
        {0, {4}}, {1, {5}}, {2, {6}}, {3, {7}}
    };

    // Inizializza la mappa degli estremi dei nuovi lati ai nuovi lati
    std::map<std::array<unsigned int, 2>, unsigned int> mappa_estremi_nuovi_lati_nuovi_lati = {
        {{1, 0}, 4}, {{0, 3}, 7}, {{3, 2}, 6}, {{2, 1}, 5}
    };

    // Inizializza i poligoni
    std::vector<unsigned int> poligoni = {0};

    unsigned int vecchio_lato_partenza = 0;
    unsigned int vecchio_lato_arrivo = 2; // Cambiato per evitare il caso degenere
    unsigned int est_taglio_partenza = 0;
    unsigned int est_taglio_arrivo = 2; // Cambiato per evitare il caso degenere
    unsigned int ID_1D = 4;

    // Chiama la funzione nuovo_poligono
    std::array<std::vector<unsigned int>, 2> result = LibraryDFN::nuovo_poligono(
        mesh, mappa_vecchi_lati_nuovi_lati, mappa_estremi_nuovi_lati_nuovi_lati,
        poligoni, 0, vecchio_lato_partenza, vecchio_lato_arrivo,
        est_taglio_partenza, est_taglio_arrivo, ID_1D
        );

    // Verifica che il risultato sia corretto
    std::vector<unsigned int> expected_vertices = {2, 0, 1};
    std::vector<unsigned int> expected_edges = {4, 4, 1};

    EXPECT_EQ(result[0], expected_vertices);
    EXPECT_EQ(result[1], expected_edges);
}

// Test per il caso degenere (triangolo)
TEST(NuovoPoligonoTest, CasoDegenere) {
    // Inizializza PolygonalMesh
    PolygonalLibrary::PolygonalMesh mesh;
    mesh.NumberCell0D = 4;
    mesh.Cell0DId = {0, 1, 2, 3};
    mesh.Cell0DCoordinates = {Eigen::Vector3d(0,0,0), Eigen::Vector3d(1,0,0), Eigen::Vector3d(1,1,0), Eigen::Vector3d(0,1,0)};
    mesh.Cell0DMarkers = {0, 0, 0, 0};

    mesh.NumberCell1D = 4;
    mesh.Cell1DId = {0, 1, 2, 3};
    mesh.Cell1DVertices = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    mesh.Cell1DMarkers = {0, 0, 0, 0};

    mesh.NumberCell2D = 1;
    mesh.Cell2DId = {0};
    mesh.Cell2DVertices = {{0, 1, 2, 3}};
    mesh.Cell2DEdges = {{0, 1, 2, 3}};

    // Inizializza la mappa dei vecchi lati ai nuovi lati
    std::map<unsigned int, std::vector<unsigned int>> mappa_vecchi_lati_nuovi_lati = {
        {0, {4}}, {1, {5}}, {2, {6}}, {3, {7}}
    };

    // Inizializza la mappa degli estremi dei nuovi lati ai nuovi lati
    std::map<std::array<unsigned int, 2>, unsigned int> mappa_estremi_nuovi_lati_nuovi_lati = {
        {{1, 0}, 4}, {{0, 3}, 7}, {{3, 2}, 6}, {{2, 1}, 5}
    };

    // Inizializza i poligoni
    std::vector<unsigned int> poligoni = {0};

    unsigned int vecchio_lato_partenza = 0;
    unsigned int vecchio_lato_arrivo = 2;
    unsigned int est_taglio_partenza = 0;
    unsigned int est_taglio_arrivo = 2;
    unsigned int ID_1D = 5;

    // Chiama la funzione nuovo_poligono
    std::array<std::vector<unsigned int>, 2> result =  LibraryDFN::nuovo_poligono(
        mesh, mappa_vecchi_lati_nuovi_lati, mappa_estremi_nuovi_lati_nuovi_lati,
        poligoni, 0, vecchio_lato_partenza, vecchio_lato_arrivo,
        est_taglio_partenza, est_taglio_arrivo, ID_1D
        );

    // Poiché il taglio passa esattamente sui vertici, il risultato atteso è che i nuovi vertici e lati siano specifici
    std::vector<unsigned int> expected_vertices = {2, 0, 1};
    std::vector<unsigned int> expected_edges = {5, 4, 1};

    EXPECT_EQ(result[0], expected_vertices);
    EXPECT_EQ(result[1], expected_edges);
}












