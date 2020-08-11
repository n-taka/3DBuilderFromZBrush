
#include <vector>
#include <string>
#include "readGoZFile.h"
#include "repairAs3DBuilder.h"
#include "writeGoZFile.h"

void main()
{
    std::string meshName;
    std::vector<std::vector<double>> V, VC;
    std::vector<double> M;
    std::vector<int> G;
    std::vector<std::vector<int>> F;
    std::vector<std::vector<std::pair<double, double>>> UV;

    FromZ::readGoZFile("./example/Dog.GoZ", meshName, V, F, UV, VC, M, G);

    std::vector<std::vector<double>> V_repair, VC_repair;
    std::vector<double> M_repair;
    std::vector<int> G_repair;
    std::vector<std::vector<int>> F_repair;
    std::vector<std::vector<std::pair<double, double>>> UV_repair;

    repairAs3DBuilder(V, F, UV, VC, M, G, V_repair, F_repair, UV_repair, VC_repair, M_repair, G_repair);

    FromZ::writeGoZFile("./example/DogNotRepairOut.GoZ", meshName, V, F, UV, VC, M, G);
    FromZ::writeGoZFile("./example/DogRepairOut.GoZ", meshName, V_repair, F_repair, UV_repair, VC_repair, M_repair, G_repair);


    return;
}