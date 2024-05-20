#pragma once

#include "StructDFN.hpp"

namespace LibraryDFN
{
    void readDFNFromFile(const std::string& filename, DFN& dfn);
    void printTraces(const DFN& dfn, const std::string& filename);
    void printTracesByFracture(const DFN& dfn, const std::string& filename);
}


