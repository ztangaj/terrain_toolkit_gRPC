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

using grpc::Channel;
using grpc::ClientContext;
using grpc::Status;
using geodesic_gRPC::HelloRequest;
using geodesic_gRPC::HelloReply;
using geodesic_gRPC::FindPathByVertexCordRequest;
using geodesic_gRPC::Path;
using geodesic_gRPC::Geodesic;

class GeodesicClient {
 public:
  GeodesicClient(std::shared_ptr<Channel> channel)
      : stub_(Geodesic::NewStub(channel)) {}

  // Assembles the client's payload, sends it and presents the response back
  // from the server.
  std::string SayHello(const std::string& user) {
    // Data we are sending to the server.
    HelloRequest request;
    request.set_name(user);

    // Container for the data we expect from the server.
    HelloReply reply;

    // Context for the client. It could be used to convey extra information to
    // the server and/or tweak certain RPC behaviors.
    ClientContext context;

    // The actual RPC.
    Status status = stub_->SayHello(&context, request, &reply);

    // Act upon its status.
    if (status.ok()) {
      return reply.message();
    } else {
      std::cout << status.error_code() << ": " << status.error_message()
                << std::endl;
      return "RPC failed";
    }
  }

  std::string FindPathByVertexCord(const std::string& user) {
    // Follows the same pattern as SayHello.
    FindPathByVertexCordRequest request;
    request.set_algo_type(user);
    Path path;
    ClientContext context;

    // Here we can use the stub's newly available method we just added.
    Status status = stub_->FindPathByVertexCord(&context, request, &path);
    if (status.ok()) {
      std::cout.precision(12);
      std::cout << path.path().Get(1);
      return path.message();
    } else {
      std::cout << status.error_code() << ": " << status.error_message()
                << std::endl;
      return "RPC failed";
    }
  }
 private:
  std::unique_ptr<Geodesic::Stub> stub_;
};

int main(int argc, char** argv) {
  // Instantiate the client. It requires a channel, out of which the actual RPCs
  // are created. This channel models a connection to an endpoint (in this case,
  // localhost at port 50051). We indicate that the channel isn't authenticated
  // (use of InsecureChannelCredentials()).
  GeodesicClient client(grpc::CreateChannel(
      "localhost:50051", grpc::InsecureChannelCredentials()));
  std::string user("world");
  std::string reply = client.SayHello(user);
  std::cout << "Greeter received: " << reply << std::endl;

  reply = client.FindPathByVertexCord(user);
  std::cout << "Greeter received: " << reply << std::endl;

  return 0;
}
