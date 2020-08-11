#include <vector>
#include <string>
#include <filesystem>
#include "readGoZFile.h"
#include "repairAs3DBuilder.h"
#include "writeGoZFile.h"

void main(int argc, char* argv[])
{
    if (argc >= 2)
    {
        std::string inputGoZFileName;
        std::string meshName;
        std::vector<std::vector<double>> V, VC;
        std::vector<double> M;
        std::vector<int> G;
        std::vector<std::vector<int>> F;
        std::vector<std::vector<std::pair<double, double>>> UV;

        inputGoZFileName = std::string(argv[1]);

        FromZ::readGoZFile(inputGoZFileName, meshName, V, F, UV, VC, M, G);

        std::string outputGoZFileName;
        std::vector<std::vector<double>> V_repair, VC_repair;
        std::vector<double> M_repair;
        std::vector<int> G_repair;
        std::vector<std::vector<int>> F_repair;
        std::vector<std::vector<std::pair<double, double>>> UV_repair;

        repairAs3DBuilder(V, F, UV, VC, M, G, V_repair, F_repair, UV_repair, VC_repair, M_repair, G_repair);

        std::filesystem::path input(inputGoZFileName);
        std::filesystem::path output = input.parent_path();
        output /= input.stem();
        output += "_repaired.GoZ";
        outputGoZFileName = output.string();

        FromZ::writeGoZFile(outputGoZFileName, meshName, V_repair, F_repair, UV_repair, VC_repair, M_repair, G_repair);
    }
    else
    {
        std::cerr << "Please drag-and-drop GoZ file to this .exe file." << std::endl;
    }

    return;
}