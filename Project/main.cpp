#include "StructDFN.hpp"
#include "Utils.hpp"
#include "Eigen/Eigen"

using namespace LibraryDFN;
using namespace Eigen;

int main()
{
    DFN mannaggia_dio;
    mannaggia_dio.numFratture=4;
    mannaggia_dio.idFratture={0,1,2,3};
    mannaggia_dio.numVertici={3,3,3,3};
    mannaggia_dio.vertici={{Vector3d(1.,1.,0.),Vector3d(1.,4.,0.),Vector3d(4.,1.,0.)},{Vector3d(0.,1.,12.),Vector3d(0.,14.,2.),Vector3d(0.,34.,1.)},{Vector3d(1.,1.,0.2),Vector3d(1.,4.,0.2),Vector3d(4.,1.,0.2)},{Vector3d(1.,1.,100.),Vector3d(1.,4.,100.),Vector3d(4.,1.,100.)}};
    std::vector<std::array<unsigned int,2>> nonloso;
    nonloso=scarta_fratture(mannaggia_dio);
    return 0;
}
