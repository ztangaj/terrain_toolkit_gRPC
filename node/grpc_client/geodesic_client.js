var PROTO_PATH = __dirname + '/../../grpc/protos/geodesic.proto';

var grpc = require('grpc');
var express = require('express')
var cors = require('cors')
var bodyParser = require('body-parser');
var multer = require('multer');

var app = express()
var forms = multer();

app.use(cors())

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

// parse application/json
app.use(bodyParser.json())
// parse application/x-www-form-urlencoded
app.use(bodyParser.urlencoded({ extended: false }))
app.use(forms.array()); 


// TODO: handle error
function getPath(query, successCallback, errorCallback) {
  var client = new geodesic_proto.Geodesic('localhost:50051', grpc.credentials.createInsecure());
  client.FindPathByVertexCord(query, function(err, response){
    successCallback(response.path);
    errorCallback(err);
  });
}


// function getTerrainInfo(query, successCallback, errorCallback){
//   var client = new geodesic_proto.Geodesic('localhost:50051', grpc.credentials.createInsecure());
//   client.GetTerrainInfo(query, function(err, response){
//     successCallback(response.info);
//     errorCallback(err);
//   });
// }

function loadModel(query, successCallback, errorCallback){
  var client = new geodesic_proto.Geodesic('localhost:50051', grpc.credentials.createInsecure());
  client.LoadModel(query, function(err, response){
    successCallback(response.info);
    errorCallback(err);
  });
}

function simplifyTerrain(query, successCallback, errorCallback){
  var client = new geodesic_proto.Geodesic('localhost:50051', grpc.credentials.createInsecure());
  client.SimplifyTerrain(query, function(err, response){
    successCallback(response.model_path);
    errorCallback(err);
  });
}

function uploadModelContent(query, successCallback, errorCallback){
  var client = new geodesic_proto.Geodesic('localhost:50051', grpc.credentials.createInsecure());
  client.UploadModelContent(query, function(err, response){
    successCallback(response.model_path);
    errorCallback(err);
  });
}

// wraps the above function call into a Promise
// and handles the callbacks with resolve and reject


function getPathAsync(query) {
  return new Promise((resolve, reject) => {
    getPath(query,(successResponse) => {
          resolve(successResponse);
      }, (errorResponse) => {
          console.log(errorResponse);
          reject(errorResponse)
      });
  });
}

// function getTerrainInfoAsync(query){
//   return new Promise((resolve, reject) => {
//     getTerrainInfo(query,(successResponse) => {
//           resolve(successResponse);
//       }, (errorResponse) => {
//           console.log(errorResponse);
//           reject(errorResponse)
//       });
//   });
// }

function loadModelAsync(query){
  return new Promise((resolve, reject) => {
    loadModel(query,(successResponse) => {
          resolve(successResponse);
      }, (errorResponse) => {
          console.log(errorResponse);
          reject(errorResponse)
      });
  });
}

function simplifyTerrainAsync(query){
  return new Promise((resolve, reject) => {
    simplifyTerrain(query,(successResponse) => {
          resolve(successResponse);
      }, (errorResponse) => {
          console.log(errorResponse);
          reject(errorResponse)
      });
  });
}

function uploadModelContentAsync(query){
  return new Promise((resolve, reject) => {
    uploadModelContent(query,(successResponse) => {
          resolve(successResponse);
      }, (errorResponse) => {
          console.log(errorResponse);
          reject(errorResponse)
      });
  });
}

// API
app.post('/getpath', function(request, response) {
  console.log("post request for path received");
  // console.log(request.body.v1);
  // console.log(request.body.v2);
  getPathAsync(request.body).then((result)=>{
    response.json({path: result.toString()});
    response.status(200).end();
  })
})

// app.get('/get-terrain-info', function(request, response) {
//   console.log("receive get request for terrain info received");
//   getTerrainInfoAsync(request.body).then((result)=>{
//     console.log(result);
//     response.json({info: result.toString()});
//     response.status(200).end();
//   })
// })

app.post('/load-model', function(request, response) {
  console.log("post request for load model received");
  loadModelAsync(request.body).then((result)=>{
    console.log(result);
    response.json({info: result.toString()});
    response.status(200).end();
  })
})

app.post('/upload-model-content', function(request, response) {
  console.log("post request for upload model received");
  uploadModelContentAsync(request.body).then((result)=>{
    console.log(result);
    response.json({model_path: result.toString()});
    response.status(200).end();
  })
})

app.post('/simplify-terrain', function(request, response) {
  console.log("post request for terrain simplification received");
  simplifyTerrainAsync(request.body).then((result)=>{
    console.log(result);
    response.json({model_path: result.toString()});
    response.status(200).end();
  })
})


port = 8080
app.listen(port)
console.log(`Listening at http://localhost:${port}`)