#include "StructDFN.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include <list>
#include <tuple>

namespace LibraryDFN{

std::vector<std::array<unsigned int,3>> triangola_frattura(DFN& disc_frac_net,unsigned int idr);
std::vector<std::array<unsigned int,2>> scarta_fratture(DFN& disc_frac_net);
//void memorizza_tracce(DFN& disc_frac_net);
}
