// cross platflam defining
#ifndef defined ( WIN32 )
#define __func__ __FUNCTION__
#endif
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include <fstream>
#include <string>
#include "critical.h"
#include <tuple>
#include <vector>

typedef std::vector<std::tuple<int,int> > Edge_list;

using namespace simplify;
// WeightMap weightmap = get(boost::edge_weight, g);  

bool exist(std::tuple<int, int> value, std::vector<std::tuple<int,int>> v){
    // return std::find(v.begin(), v.end(),value)!=v.end();
    for(auto it = v.begin(); it != v.end(); ++it) {
        if(std::get<0>(value) == std::get<0>(*it) && std::get<1>(value) == std::get<1>(*it)){
            return true;
        }
        if(std::get<0>(value) == std::get<1>(*it) && std::get<1>(value) == std::get<0>(*it)){
            return true;
        }
    }
    return false;
}

double distance(Point3 v1, Point3 v2){
    double d[] = {abs(v1.x-v2.x), abs(v1.y-v2.y), abs(v1.z-v2.z)};
    // if (d[0] < d[1]) std::swap(d[0],d[1]);
    // if (d[0] < d[2]) std::swap(d[0],d[2]);
    // double distance = d[0] * sqrt(1.0 + d[1]/d[0] + d[2]/d[0]);
    return sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
}

// we assume there is not comment in OFF file
bool generate_graph(std::string in_path, std::string out_path){
    // std::ifstream infile("/etc/terrain_toolkit/models/off/small_terrain.off");
    // std::ofstream outfile("/etc/terrain_toolkit/models/graph/small_terrain.graph");

    std::ifstream infile(in_path);
    std::ofstream outfile(out_path);

    Edge_list edge_list;

    if (infile.is_open() && outfile.is_open()) {
        std::string line;
        int num_v, num_f, num_e;
        double x,y,z;
        int k,v1,v2,v3;
        infile >> line;
        if(line.compare("OFF")!=0){
            std::cout<<"not off file!"<<std::endl;
        }
        infile >> num_v >> num_f >> num_e;
        outfile << "g "<< num_v << " "<< num_e << std::endl;
        // read vertex
        for(int i=0; i<num_v; i++){
            infile >> x >> y >> z;
            vertex_props.push_back( Point3{ x,y,z } );
            outfile << i << " " << x << " " << y << " " << z << std::endl;
        }
        // read face
        for(int i=0;i<num_f;i++){
            infile >> k >> v1 >> v2 >> v3;

            if(!exist(std::make_tuple(v1,v2), edge_list)){
                edge_list.push_back(std::make_tuple(v1,v2));
                outfile << v1 << " " << v2 << " " << distance(vertex_props[v1],vertex_props[v2]) << std::endl;
            }

            if(!exist(std::make_tuple(v1,v3), edge_list)){
                edge_list.push_back(std::make_tuple(v1,v3));
                outfile << v1 << " " << v3 << " " << distance(vertex_props[v1],vertex_props[v3]) << std::endl;
            }

            if(!exist(std::make_tuple(v2,v3), edge_list)){
                edge_list.push_back(std::make_tuple(v2,v3));
                outfile << v2 << " " << v3 << " " << distance(vertex_props[v2],vertex_props[v3]) << std::endl;
            }
        }
        infile.close();
        outfile.close();
        return true;
    }
    else{
        std::cout<<"file to open file" <<std::endl;
        return false;
    }
}

