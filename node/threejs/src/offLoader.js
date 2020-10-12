class OffLoader {
    /*  parser rules
        1. if see #, skip to next line
        2. first line, OFF
        3. second line, num_vertices num_faces num_edges
        4. vertices line, x y z
        5. faces line, dim v1 v2 v3
        TODO:
        we only consider case dim==3
        we assume no more comment after the num info line
        // reference https://en.wikipedia.org/wiki/OFF_(file_format)
    */

    // this method will return model =  {"v":vertices, "f":faces, "z":{"min":minz, "max":maxz}};
    static loadModel(content){
        var vertices = [];
        var faces = [];
        const allLines =  content.split(/\r\n|\n/);
        var validLineCount = 0;
        var vertexCount, faceCount, edgeCount;
        
        // Reading line by line
        var vertexLineStart = 0;
        for(vertexLineStart; vertexLineStart<allLines.length; vertexLineStart++){
            var line = allLines[vertexLineStart];
            if(line.startsWith("#")){  
                continue;
            }           
            if(validLineCount==0){
                // OFF
                if(!(line == "OFF")){
                    console.log("Wrong data format: " + line);
                    break;                
                }
                validLineCount += 1;
            }
            else if(validLineCount==1){
                // count info
                var infos = line.split(" ");
                vertexCount = parseInt(infos[0], 10);
                faceCount = parseInt(infos[1], 10);
                edgeCount = parseInt(infos[2], 10);
                validLineCount += 1;
            }               
            else if(validLineCount==2){
                break;
            }
        }
        // using more loops, we can get rid of if branches
        // get list of vertices
        var x,y,z;
        var maxz = -1e9;
        var minz = 1e9;
        for(var i=0; i<vertexCount; i++){
            var line = allLines[i+vertexLineStart];
            var infos = line.split(" ");
            x = parseFloat(infos[0], 10);
            y = parseFloat(infos[1], 10);
            z = parseFloat(infos[2], 10);
            if(z>maxz){
                maxz=z;
            }
            if(z<minz){
                minz=z;
            }
            var v = new THREE.Vector3(x, y, z);
            vertices.push(v);
        }
        // get list of faces
        var v1, v2, v3;
        var faceLineStart = vertexLineStart + vertexCount;
    
        for(var i=0; i<faceCount; i++){
            var line = allLines[i+faceLineStart];
            var infos = line.split(" ");
            v1 = parseInt(infos[1], 10);
            v2 = parseInt(infos[2], 10);
            v3 = parseInt(infos[3], 10);        
            faces.push([v1,v2,v3]);
            // edges implementation
            // var e1 = new VertexPair(v1,v2);
            // var e2 = new VertexPair(v1,v3);
            // var e3 = new VertexPair(v2,v3);
            // var temp = [e1,e2,e3];
            // temp.forEach(element => {
            //     if(!edge_exist(edges, element)){
            //         edges.push(element);
            //     }
            // });    
        }
        var model =  {"v":vertices, "f":faces, "z":{"min":minz, "max":maxz}};
        console.log(model);
        return model;
    }

    static visualizeModel(camera, controls, scene, model, material){
        var vertices = model.v;
        var faces= model.f;
        var group = new THREE.Group();
    
        faces.forEach(f => {
            var geometry = new THREE.BufferGeometry().setFromPoints([vertices[f[0]], vertices[f[1]] , vertices[f[2]], vertices[f[0]]]);
            var line = new THREE.Line( geometry, material );     
            group.add(line);
        });
        var boundingbox = new THREE.Box3();
        boundingbox.setFromObject(group);
        var center = boundingbox.getCenter();
        var sceneRadiusForCamera = Math.max(
            boundingbox.max.y - boundingbox.min.y,
            boundingbox.max.z - boundingbox.min.z,
            boundingbox.max.x - boundingbox.min.x
        )/2 * (1 + Math.sqrt(5));
            
        camera.position.set(center.x, center.y, center.z + sceneRadiusForCamera);
        controls.target.set(center.x, center.y, center.z);
        scene.add(group);
        controls.update();
    }

    static exportModel(vertices, faces){
        var res = "";
        res += "OFF\n";
        // TODO: in some format we dont have edge count
        res += vertices.length + " " + faces.length + " " + faces.length + "\n";
        for(var i=0; i<vertices.length; i++){
            var v = vertices[i];
            res += v.x + " " + v.y + " " + v.z + "\n";
        }
        for(var i=0; i<faces.length; i++){
            var f = faces[i];
            res += "3 " + f[0] + " " + f[1] + " " + f[2] + "\n";
        }
        return res;
    }
}