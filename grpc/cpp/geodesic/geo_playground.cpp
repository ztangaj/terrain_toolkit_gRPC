#include "distance.h"
#include<sstream>
#include<unistd.h>
#include <sys/resource.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/times.h>
#include "geodesic_algorithm_subdivision.h"
#include "geodesic_mesh_elements.h"
#include "geodesic_memory.h"
#include "geodesic_constants_and_simple_functions.h"

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

    geodesic::GeodesicAlgorithmExact algorithm(&mesh);	//create exact algorithm for the mesh

    // geodesic::SurfacePoint source(&mesh.vertices()[atol(argv[2])]);

    // point, mesh, find nearest vertex
    // geodesic::MeshElementBase base = geodesic::MeshElementBase();
    // geodesic::SurfacePoint source(&base, 589230,5213190,1406, geodesic::VERTEX);
    for(int i=0; i<mesh.vertices().size();i++){
        geodesic::Vertex v = mesh.vertices()[i];
        if(v.x()==589230 && v.y()==5213190){
            std::cout<<"found match vertex"<<std::endl;
        }
    }
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
