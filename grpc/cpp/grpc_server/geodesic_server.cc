/*
 *
 * Copyright 2015 gRPC authors.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <iostream>
#include <memory>
#include <string>

#include <grpcpp/grpcpp.h>
#include "geodesic.grpc.pb.h"

#include<sstream>
#include<unistd.h>
#include <sys/resource.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <sys/times.h>

#include "distance.h"
#include "geodesic_algorithm_subdivision.h"

using grpc::Server;
using grpc::ServerBuilder;
using grpc::ServerContext;
using grpc::Status;
using geodesic_gRPC::HelloRequest;
using geodesic_gRPC::HelloReply;
using geodesic_gRPC::FindPathByVertexCordRequest;
using geodesic_gRPC::Path;
using geodesic_gRPC::Geodesic;

geodesic::SurfacePoint findVertexByCordStr(std::string cordStr, geodesic::Mesh* mesh){
  std::string delimiter = ",";
  size_t pos = 0;
  std::string token;
  std::vector<_Float32> cords;
  geodesic::SurfacePoint result;

  while ((pos = cordStr.find(delimiter)) != std::string::npos) {
      token = cordStr.substr(0, pos);
      std::cout << token << std::endl;
      cordStr.erase(0, pos + delimiter.length());
      cords.push_back(std::stof(token));
  }

  for(int i=0; i<mesh->vertices().size();i++){
      geodesic::SurfacePoint v(&mesh->vertices()[i]);
      if(v.x()==cords[0] && v.y()==cords[1]){
          std::cout<<"found match vertex"<<std::endl;        
          result = v;
      }
  }
  return result;
}

// Logic and data behind the server's behavior.
class GeodesicServiceImpl final : public Geodesic::Service {
  Status SayHello(ServerContext* context, const HelloRequest* request,
                  HelloReply* reply) override {
    std::string prefix("Hello ");
    reply->set_message(prefix + request->name());
    return Status::OK;
  }
  Status FindPathByVertexCord(ServerContext* context, const FindPathByVertexCordRequest* request,
                       Path* reply) override {
    bool success = geodesic::read_mesh_from_file("small_terrain.off",points,faces);

    if(!success)
    {
      reply->set_message("something is wrong with the input file");
      return Status::CANCELLED;
    }

    mesh.initialize_mesh_data(points, faces);		//create internal mesh data structure including edges

    geodesic::GeodesicAlgorithmExact algorithm(&mesh);	//create exact algorithm for the mesh

    // std::cout << request->algo_type() << std::endl;
    // std::cout << request->v1() << std::endl;
    // std::cout << request->v2() << std::endl;
    geodesic::SurfacePoint source = findVertexByCordStr(request->v1(), &mesh);
    geodesic::SurfacePoint destination = findVertexByCordStr(request->v2(), &mesh);

    // geodesic::SurfacePoint source(&mesh.vertices()[5]);
    // geodesic::SurfacePoint destination(&mesh.vertices()[500]);

    std::vector<geodesic::SurfacePoint> path;
    algorithm.geodesic(source, destination, path);

    for(unsigned i = 0; i<path.size(); ++i)
    {
        geodesic::SurfacePoint& s = path[i];
        
        std::cout << s.x() << "\t" << s.y() << "\t" << s.z() << std::endl;
        reply->add_path(s.x());
        reply->add_path(s.y());
        reply->add_path(s.z());
    }
    reply->set_message("success");
    return Status::OK;
  }
};

void RunServer() {
  std::string server_address("0.0.0.0:50051");
  GeodesicServiceImpl service;

  ServerBuilder builder;
  // Listen on the given address without any authentication mechanism.
  builder.AddListeningPort(server_address, grpc::InsecureServerCredentials());
  // Register "service" as the instance through which we'll communicate with
  // clients. In this case it corresponds to an *synchronous* service.
  builder.RegisterService(&service);
  // Finally assemble the server.
  std::unique_ptr<Server> server(builder.BuildAndStart());
  std::cout << "Server listening on " << server_address << std::endl;

  // Wait for the server to shutdown. Note that some other thread must be
  // responsible for shutting down the server for this call to ever return.
  server->Wait();
}

int main(int argc, char** argv) {
  RunServer();

  return 0;
}
