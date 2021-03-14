./grpc/cpp/grpc_server/cmake/build/geodesic_server >> grpc_server.log 2>&1 &
node node/grpc_client/geodesic_client.js >> node_server.log 2>&1 &
