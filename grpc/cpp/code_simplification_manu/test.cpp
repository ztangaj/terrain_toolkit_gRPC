#include "critical.h"
using namespace simplify;

int main(int argc, char* argv[])
{
    std::string offpath = "//home/tongcheng/Documents/terrain_toolkit_web/node/threejs/src/models/off/small_terrain.off.graph";
    float beta = 1.8;
    std::string outpath = "test4.off";
    generateOff(offpath, beta, outpath);
}