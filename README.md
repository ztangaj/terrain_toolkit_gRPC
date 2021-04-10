# terrain_toolkit_gRPC Tutorial
Based on Ubuntu 20.04

# Pre-requisites

## Some libraries

```ruby
 $ sudo apt-get install build-essential autoconf libtool pkg-config
 $ sudo apt-get install cmake
```

## Clone from github
Clone this github repository and the grpc sbumodule
```ruby
 $ git clone https://github.com/ztangaj/terrain_toolkit_gRPC.git
 $ cd terrain_toolkit_gRPC
 $ git submodule update --init
```
You can see that there are two folders _/grpc_ and _/node_. 
_/grpc_ is actually an terrain algorithm written by c++, we should create a gRPC server here. And the _/node_ is an node.js project with gRPC client. 

## Change Path
in vs code, search
```
home/tongcheng/Documents/terrain_toolkit_web
```
replace to your own path (eg. 'xxx/terrain_toolkit_gRPC')

# GPRC Server

## Clone the repository in _grpc_
```ruby
 $ cd grpc
 $ git clone https://github.com/grpc/grpc
 $ cd grpc
 $ git submodule update --init
```
## Build from source
```ruby
 $ mkdir -p cmake/build
 $ cd cmake/build
 $ cmake ../..
 $ make
```

## Install after build
*If you have any problem, see the Troubleshooting*
```ruby
# NOTE: all of gRPC's dependencies need to be already installed
$ cd ../../../cpp/grpc_server/
$ cd cmake/build
$ cmake ../.. -DgRPC_INSTALL=ON                \
              -DCMAKE_BUILD_TYPE=Release       \
              -DgRPC_ABSL_PROVIDER=package     \
              -DgRPC_CARES_PROVIDER=package    \
              -DgRPC_PROTOBUF_PROVIDER=package \
              -DgRPC_RE2_PROVIDER=package      \
              -DgRPC_SSL_PROVIDER=package      \
              -DgRPC_ZLIB_PROVIDER=package
$ make
```

## Install boost-1.65.1
Since apt-get do not provide old version of boost, we should install it from SourceForge
Click the following link to install in your own local path (eg. _/path/to_)
https://www.boost.org/doc/libs/1_65_1/more/getting_started/unix-variants.html
```ruby
 $ tar --bzip2 -xf /path/to/boost_1_65_1.tar.bz2
 $ cd /path/to/boost_1_65_1
 $ ./bootstrap.sh
 $ ./b2
```
## Set the path of boost-1.65.1
(change _/path/to_ to your own path that contains boost_1_65_1, for me, it is _/home/tongcheng/Downloads/_)
```ruby
 $ cd /terrain_toolkit_gRPC/grpc/cpp/grpc_server/cmake/build
 $ unset LD_LIBRARY_PATH
 $ export LD_LIBRARY_PATH=/path/to/boost_1_65_1/stage/lib:$LD_LIBRARY_PATHD
```

## Check ldd to make sure all of the required libraries found
```ruby
 $ ldd geodesic_server
```

## Start the gRPC server!
```ruby
 $ ./geodesic_server
```


## Troubleshooting
### Can't use protobuf in cmakelists.txt
https://stackoverflow.com/questions/41573702/cant-use-protobuf-in-cmakelists-txt
error: 
```ruby
CMake Error at CMakeLists.txt:9 (find_package):
  Could not find a package configuration file provided by "protobuf" with any
  of the following names:

    protobufConfig.cmake
    protobuf-config.cmake

  Add the installation prefix of "protobuf" to CMAKE_PREFIX_PATH or set
  "protobuf_DIR" to a directory containing one of the above files.  If
  "protobuf" provides a separate development package or SDK, be sure it has
  been installed.
```

Find: 
```ruby
find_package(protobuf CONFIG REQUIRED)
```
Replace to: 
```ruby
find_package(Protobuf REQUIRED)
```


## Reference: 
https://chromium.googlesource.com/external/github.com/grpc/grpc/+/HEAD/BUILDING.md
https://grpc.io/docs/languages/cpp/quickstart/
https://gitee.com/lancelotghx/terrain_toolkit_web/tree/0a186d15f4f963dbb8797c0c9d8f8d0520a509e1/node
https://github.com/grpc/grpc/issues/11378

# GPRC Client
*You should start the client in a new terminal*

## Install some required npm packages
```ruby
 $ cd node
 $ npm i
 $ sudo npm install grpc
 $ cd grpc client
 $ node ./geodesic_client.js
 $ npm install express
 $ npm install @grpc/proto-loader
```

## Start the client
```ruby
 $ node ./geodesic client.js
```

# HTTP Server
First, back to the terrain_toolkit_gRPC folder
```ruby
 $ sudo npm install --global http-server
 $ http-server ./ 
```

# Start
Goto http://127.0.0.1:8080/node/threejs/src/terrain_toolkit.html