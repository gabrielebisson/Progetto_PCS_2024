#pragma once

#include <gtest/gtest.h>
#include <fstream>
#include "StructDFN.hpp"
#include "Utils.hpp"
#include <Eigen/Eigen>

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

    std::vector<Eigen::Vector3d> versori_normali = {Eigen::Vector3d(0, 0, 1), Eigen::Vector3d(0, 0, 1)};
    double tol = 0.1;

    auto result = scarta_fratture(dfn, versori_normali, tol);
    EXPECT_TRUE(result.empty());
}

// Test per il caso in cui le fratture si intersecano
TEST(ScartaFrattureTest, Intersection) {
    LibraryDFN::DFN dfn;
    dfn.numFratture = 2;
    dfn.vertici = {
        {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(2, 0, 0), Eigen::Vector3d(2, 2, 0), Eigen::Vector3d(0, 2, 0)},
        {Eigen::Vector3d(1, 1, -1), Eigen::Vector3d(3, 1, -1), Eigen::Vector3d(3, 1, 1), Eigen::Vector3d(1, 1, 1)}
    };
    dfn.numVertici = {4, 4};

    std::vector<Eigen::Vector3d> versori_normali = {Eigen::Vector3d(0, 0, 1), Eigen::Vector3d(0, 1, 0)};
    double tol = 0.1;

    auto result = scarta_fratture(dfn, versori_normali, tol);

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

    std::vector<Eigen::Vector3d> versori_normali = {Eigen::Vector3d(0, 0, 1), Eigen::Vector3d(0, 0, 1)};
    double tol = 0.1;

    auto result = scarta_fratture(dfn, versori_normali, tol);
    EXPECT_TRUE(result.empty());
}
