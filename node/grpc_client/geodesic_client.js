var PROTO_PATH = __dirname + '/../../grpc/protos/geodesic.proto';

var fs = require('fs');
var grpc = require('grpc');
var express = require('express');
var cors = require('cors');
var bodyParser = require('body-parser');
var multer = require('multer');
var upload = multer({ dest: "upload/" })

var app = express()
// var forms = multer();

app.use(cors());

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

var SRC_DIR = "../threejs/src/";

// parse application/json
app.use(bodyParser.json())
// parse application/x-www-form-urlencoded
app.use(bodyParser.urlencoded({ extended: false }))
// app.use(forms.array()); 

// TODO: handle error
function getPath(query, successCallback, errorCallback) {
  var client = new geodesic_proto.Geodesic('localhost:50051', grpc.credentials.createInsecure());
  client.FindPathByVertexCord(query, function(err, response){
    if(response){      
      if(response.path.length>0){
        successCallback(response.path);
      }
      errorCallback(response.message);
    }
    else{
      errorCallback("Server no response");
    }
  });
}

function loadModel(query, successCallback, errorCallback){
  var client = new geodesic_proto.Geodesic('localhost:50051', grpc.credentials.createInsecure());
  client.LoadModel(query, function(err, response){
    if(response){
      successCallback(response.info);
    }
    else{
      errorCallback(err);
    }
  });
}

function simplifyTerrain(query, successCallback, errorCallback){
  var client = new geodesic_proto.Geodesic('localhost:50051', grpc.credentials.createInsecure());
  client.SimplifyTerrain(query, function(err, response){
    if(response){
      if(response.model_path != ""){
        successCallback(response.model_path);
      }
      errorCallback(response.message);
    }
    else{
      errorCallback("Server no response");
    }
  });
}

// wraps the above function call into a Promise
// and handles the callbacks with resolve and reject


function getPathAsync(query) {
  return new Promise((resolve, reject) => {
    getPath(query,(successResponse) => {
          resolve(successResponse);
      }, (errorResponse) => {
        reject(new Error(errorResponse));
      });
  });
}

function loadModelAsync(query){
  return new Promise((resolve, reject) => {
    loadModel(query,(successResponse) => {
          resolve(successResponse);
      }, (errorResponse) => {
          // console.log(errorResponse);
          reject(new Error(errorResponse));
      });
  });
}

function simplifyTerrainAsync(query){
  return new Promise((resolve, reject) => {
    simplifyTerrain(query,(successResponse) => {
          resolve(successResponse);
      }, (errorResponse) => {
          // console.log(errorResponse);
          reject(new Error(errorResponse));
      });
  });
}

// API
app.post('/getpath', function(request, response) {
  console.log("post request for find path received");
  // console.log(request.body.v1);
  // console.log(request.body.v2);
  getPathAsync(request.body)
  .then((result)=>{
    response.json({path: result.toString()});
    response.status(200).end();
  })
  .catch((err)=>{
    console.log(err);
    response.statusMessage = err;
    response.status(500).end();
  })
})

app.post('/load-model', function(request, response) {
  console.log("post request for load model received");
  loadModelAsync(request.body)
  .then((result)=>{
    console.log(result);
    response.json({info: result.toString()});
    response.status(200).end();
  })
  .catch((err)=>{
    response.status(500).end();
    console.log(err);
  })
})

app.post('/upload-model-content', upload.single("model_file"), function(request, response) {
  console.log("post request for upload model received");
  // var buf = fs.readFileSync(SRC_DIR + request.body.model_path);
  var buf = fs.readFileSync(request.file.path);
  var data = buf.toString();
  const allLines =  data.split(/\r\n|\n/);
  console.log("line count: ", allLines.length);
  var client = new geodesic_proto.Geodesic('localhost:50051', grpc.credentials.createInsecure());
  var call = client.UploadModelContent(function(error, res) {
    if (error) {
      response.status(500).end();
      console.log(error);
    }
    response.json({model_path: res.model_path});
    response.status(200).end();
  });

  for (var i = 0; i < allLines.length; i++) {
    call.write({
      content: allLines[i]
    });
  }
  call.end();
})

app.post('/simplify-terrain', multer().array(), function(request, response) {
  console.log("post request for terrain simplification received");
  simplifyTerrainAsync(request.body)
  .then((result)=>{
    console.log(result);
    response.json({model_path: result.toString()});
    response.status(200).end();
  })
  .catch((err)=>{
    console.log(err);
    response.statusMessage = err;
    response.status(500).end();
  })
})


port = 50055
app.listen(port)
console.log(`Listening at http://localhost:${port}`)