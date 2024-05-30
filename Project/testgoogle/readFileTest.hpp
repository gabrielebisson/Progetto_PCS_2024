#pragma once

#include <gtest/gtest.h>
#include "StructDFN.hpp"
#include "Utils.hpp"

// Test per la lettura del DFN da un file
TEST(ReadDFNFromFileTest, ValidInputFile) {
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

