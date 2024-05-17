#pragma once

#include "StructDFN.hpp"

namespace LibraryDFN
{
    DFN readDFNFromFile(const std::string& filename);
    void printTraces(const DFN& dfn, const std::string& filename);
    void printTracesByFracture(const DFN& dfn, const std::string& filename);
}


