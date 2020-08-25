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

// let's say this is the API function with two callbacks,
// one for success and the other for error

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

app.post('/getpath', function(request, response) {
  console.log("post request for path received");
  console.log(request.body.v1);
  console.log(request.body.v2);

  getPathAsync(request.body).then((result)=>{
    // response.writeHead(200, {'Content-Type': 'text/html'});
    console.log(result);
    response.json({path: result.toString()});
    response.status(200).end();
  })
})


port = 8080
app.listen(port)
console.log(`Listening at http://localhost:${port}`)