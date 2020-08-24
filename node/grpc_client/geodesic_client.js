var PROTO_PATH = __dirname + '/../../grpc/protos/geodesic.proto';

var grpc = require('grpc');
const http = require('http');
var express = require('express')
var fs = require('fs')
var cors = require('cors')

var app = express()
app.use(cors())
// app.use(express.bodyParser())

var protoLoader = require('@grpc/proto-loader');
const { SSL_OP_EPHEMERAL_RSA } = require('constants');
var packageDefinition = protoLoader.loadSync(
    PROTO_PATH,
    {keepCase: true,
     longs: String,
     enums: String,
     defaults: true,
     oneofs: true
    });
var geodesic_proto = grpc.loadPackageDefinition(packageDefinition).geodesic_gRPC;

// let's say this is the API function with two callbacks,
// one for success and the other for error

// TODO: handle error
function getPath(query, successCallback, errorCallback) {
  var client = new geodesic_proto.Geodesic('localhost:50051', grpc.credentials.createInsecure());
  client.FindPathByVertexID({algo_type: query}, function(err, response){
    successCallback(response.path);
  });
}

// myFunction wraps the above API call into a Promise
// and handles the callbacks with resolve and reject
function getPathAsync(query) {
  return new Promise((resolve, reject) => {
    getPath(query,(successResponse) => {
          resolve(successResponse);
      }, (errorResponse) => {
          reject(errorResponse)
      });
  });
}

// function startServer(){
//   http.createServer((request, response) => {
//     request.on('error', (err) => {
//       console.error(err);
//       response.statusCode = 400;
//       response.end();
//     });
//     response.on('error', (err) => {
//       console.error(err);
//     });
//     if (request.method === 'POST' && request.url === '/getpath') {
//       console.log("post request for path received");
//       getPathAsync("your query").then((result)=>{
//         response.setHeader("Access-Control-Allow-Origin", "*");
//         response.setHeader('Access-Control-Allow-Methods', 'GET, POST, OPTIONS, PUT, PATCH, DELETE');
//         response.setHeader("Access-Control-Allow-Headers", "*");
//         response.body = result.toString();
//         response.statusCode = 200;
//         response.write(response.body);
//         response.end();
//       });            
//     } else {
//       response.statusCode = 404;
//       response.end();
//     }
//   }).listen(8080);
// }

// function main() {
//   startServer();
// }

// main();

app.post('/getpath', function(request, response) {
  console.log("post request for path received");
  // console.dir(request.body);
  // response.writeHead(200, {'Content-Type': 'text/html'});
  // response.end('thanks');
  getPathAsync("input").then((result)=>{
    // response.writeHead(200, {'Content-Type': 'text/html'});
    console.log(result);
    response.json({path: result.toString()});
    response.status(200).end();
  })
})


port = 8080
app.listen(port)
console.log(`Listening at http://localhost:${port}`)