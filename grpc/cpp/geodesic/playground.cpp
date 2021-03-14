#include "distance.h"
#include<sstream>
#include<unistd.h>
#include <sys/resource.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/times.h>
#include <string>
#include "geodesic_algorithm_subdivision.h"
#include "geodesic_mesh_elements.h"
#include "geodesic_memory.h"
#include "geodesic_constants_and_simple_functions.h"

std::string printInfo(geodesic::Mesh mesh){
    std::string json = "";
    std::cout << "===================main===============" << std::endl;

    json += "{";
    json += "\"vertices\": " + std::to_string(mesh.vertices().size()) + ",";
    json += "\"edges\": " + std::to_string(mesh.edges().size())+ ",";
    json += "\"faces\": " + std::to_string(mesh.faces().size())+ ",";

    double minx = 1e100;
	double maxx = -1e100;
	double miny = 1e100;
	double maxy = -1e100;
	double minz = 1e100;
	double maxz = -1e100;
    for(unsigned i=0; i<mesh.vertices().size(); ++i)
	{
		geodesic::Vertex& v = mesh.vertices()[i];
		minx = std::min(minx, v.x());		
		maxx = std::max(maxx, v.x());
		miny = std::min(miny, v.y());
		maxy = std::max(maxy, v.y());
		minz = std::min(minz, v.z());
		maxz = std::max(maxz, v.z());
	}

    json += "\"minx\": " + std::to_string(minx) + ",";
    json += "\"miny\": " + std::to_string(miny) + ",";
    json += "\"minz\": " + std::to_string(minz) + ",";
    json += "\"maxx\": " + std::to_string(maxx) + ",";
    json += "\"maxy\": " + std::to_string(maxy) + ",";
    json += "\"maxz\": " + std::to_string(maxz);
    json += "}";

    std::cout << json << std::endl;
    return json;
}

int main(int argc, char **argv)
{
	if(argc < 4)
	{
		std::cout << "usage: mesh_file_name source_idx target_idx" << std::endl;
		return 0;
	}

	// bool success = geodesic::read_mesh_from_file(argv[1],points,faces);
    bool success = geodesic::read_mesh_from_file("small_terrain.off",points,faces);

	if(!success)
	{
		std::cout << "something is wrong with the input file" << std::endl;
		return 0;
	}
//    sscanf(argv[2], "%f", &s );
    s = atof(argv[2]);
	mesh.initialize_mesh_data(points, faces);		//create internal mesh data structure including edges
    // std::cout << "info: " << mesh.vertices().size()  << mesh.faces().size() << mesh.edges().size()  << std::endl;
    printInfo(mesh);
    geodesic::GeodesicAlgorithmExact algorithm(&mesh);	//create exact algorithm for the mesh

    // geodesic::SurfacePoint source(&mesh.vertices()[atol(argv[2])]);

    // point, mesh, find nearest vertex
    // geodesic::MeshElementBase base = geodesic::MeshElementBase();
    // geodesic::SurfacePoint source(&base, 589230,5213190,1406, geodesic::VERTEX);
    // for(int i=0; i<mesh.vertices().size();i++){
    //     geodesic::Vertex v = mesh.vertices()[i];
    //     if(v.x()==589230 && v.y()==5213190){
    //         std::cout<<"found match vertex"<<std::endl;
    //     }
    // }
    // std::cout << source.x() << "\t" << source.y() << "\t" << source.z() << std::endl;
    // std::vector<geodesic::Vertex*> storage;

    // mesh.closest_vertices(&source, &storage);
    // for(unsigned i = 0; i<storage.size(); ++i){
    //     geodesic::Vertex& s = *storage[i];
    //     std::cout << s.x() << "\t" << s.y() << "\t" << s.z() << std::endl;
    
    // }
    // geodesic::SurfacePoint destination(&mesh.vertices()[atol(argv[3])]);
    // std::vector<geodesic::SurfacePoint> path;

    // algorithm.geodesic(source, destination, path);

    // print_info_about_path(path);
    // std::cout.width(12);
    // std::cout.precision(10);
    // for(unsigned i = 0; i<path.size(); ++i)
    // {
    //     geodesic::SurfacePoint& s = path[i];
        
    //     std::cout << s.x() << "\t" << s.y() << "\t" << s.z() << std::endl;
    // }
}
