# terrain_toolkit_gRPC
fix some bugs in the original code and add readme file

# first, I highly recommand you to use Ubuntu 20.04

## Pre-requisites

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

## Clone the repository in _grpc_
```ruby
 $ cd grpc
 $ git clone https://github.com/grpc/grpc
 $ cd grpc
 $ git submodule update --init
 $ cd ../
```
## Build from source
```ruby
 $ mkdir -p cmake/build
 $ cd cmake/build
 $ cmake ../..
 $ make
```

## Install boost-1.65.1
Since apt-get do not provide old version of boost, we should install it from SourceForge
Click the following link to install in your own local path (eg. _path/to_)
https://www.boost.org/doc/libs/1_65_1/more/getting_started/unix-variants.html
```ruby
 $ tar --bzip2 -xf /path/to/boost_1_65_1.tar.bz2
```


## Troubleshooting
### Can't use protobuf in cmakelists.txt
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
reference: https://stackoverflow.com/questions/41573702/cant-use-protobuf-in-cmakelists-txt

reference: https://chromium.googlesource.com/external/github.com/grpc/grpc/+/HEAD/BUILDING.md
