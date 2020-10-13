class ObjLoader {
    /*  parser rules
        # List of geometric vertices, with (x, y, z [,w]) coordinates, w is optional and defaults to 1.0.
        v 0.123 0.234 0.345 1.0
        v ...
        ...
        # List of texture coordinates, in (u, [,v ,w]) coordinates, these will vary between 0 and 1. v, w are optional and default to 0.
        vt 0.500 1 [0]
        vt ...
        ...
        # List of vertex normals in (x,y,z) form; normals might not be unit vectors.
        vn 0.707 0.000 0.707
        vn ...
        ...
        # Parameter space vertices in ( u [,v] [,w] ) form; free form geometry statement ( see below )
        vp 0.310000 3.210000 2.100000
        vp ...
        ...
        # Polygonal face element (v/vt/vn)
        f 1 2 3
        f 3/1 4/2 5/3
        f 6/4/1 3/5/3 7/6/5
        f 7//1 8//2 9//3
        f ...
        ...
        # Line element (see below)
        l 5 8 1 2 4 9

        https://en.wikipedia.org/wiki/Wavefront_.obj_file
    */

    // this method will return model =  {"v":vertices, "f":faces, "z":{"min":minz, "max":maxz}};
    static loadModel(content){
        var vertices = [];
        var faces = [];
        const allLines =  content.split(/\r\n|\n/);
        
        // Reading line by line
        var x,y,z;
        var maxz = -1e9;
        var minz = 1e9;
        var v1, v2, v3;

        for(var i=0; i<allLines.length; i++){    
            var line = allLines[i];        
            if(line.startsWith("v")){  
                // parse vertex
                var infos = line.split(" ");
                x = parseFloat(infos[1], 10);
                y = parseFloat(infos[2], 10);
                z = parseFloat(infos[3], 10);
                if(z>maxz){
                    maxz=z;
                }
                if(z<minz){
                    minz=z;
                }
                var v = new THREE.Vector3(x, y, z);
                vertices.push(v);
            }   
            if(line.startsWith("f")){  
                // parse face
                var infos = line.split(" ");
                v1 = parseInt(infos[1].split("/")[0], 10)-1;
                v2 = parseInt(infos[2].split("/")[0], 10)-1;
                v3 = parseInt(infos[3].split("/")[0], 10)-1;
                faces.push([v1,v2,v3]); 
            }           
        }
  
        var model =  {"v":vertices, "f":faces, "z":{"min":minz, "max":maxz}};
        console.log(model);
        return model;
    }    

    static exportModel(vertices, faces){
        var res = "# Model exported by WebTerrain\n";
        for(var i=0; i<vertices.length; i++){
            var v = vertices[i];
            res += "v " + v.x + " " + v.y + " " + v.z + "\n";
        }
        for(var i=0; i<faces.length; i++){
            var f = faces[i];
            res += "f ".concat(f[0]+1, " ", f[1]+1, " ", f[2]+1, "\n");
        }
        return res;
    }
}