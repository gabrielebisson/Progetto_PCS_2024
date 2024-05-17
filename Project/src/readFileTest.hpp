#include <gtest/gtest.h>
#include "Utils.hpp"
#include "StructDFN.hpp"

// Test per verificare il comportamento quando il file non esiste
TEST(ReadDFNTest, NonExistentFile)
{
    std::string filename = "DFN/non_existent_file.txt";
    LibraryDFN::DFN dfn = LibraryDFN::readDFNFromFile(filename);
    // Verifica che il DFN restituito sia vuoto
    EXPECT_EQ(dfn.numFratture, 0);
}

// Test per verificare il comportamento quando il file esiste e contiene dati validi
TEST(ReadDFNTest, ValidFile)
{
    std::string filename = "DFN/correctFile";
    // Assumi che valid_file.txt contenga dati di test validi
    LibraryDFN::DFN dfn = LibraryDFN::readDFNFromFile(filename);
    // Verifica che il numero di fratture lette sia corretto
    EXPECT_EQ(dfn.numFratture, 10);
    // Puoi aggiungere ulteriori asserzioni per verificare la correttezza dei dati letti
}

