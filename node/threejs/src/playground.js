class PlyLoader {
    /*  
        The header of both ASCII and binary files is ASCII text. Only the numerical data that follows the header is different between the two versions. The header always starts with a "magic number", a line containing
        ply
        which identifies the file as a PLY file. The second line indicates which variation of the PLY format this is. It should be one of:
        format ascii 1.0
        format binary_little_endian 1.0
        format binary_big_endian 1.0
        Future versions of the standard will change the revision number at the end - but 1.0 is the only version currently in use.
        Comments may be placed in the header by using the word comment at the start of the line. Everything from there until the end of the line should then be ignored. e.g.:
        comment This is a comment!
        The 'element' keyword introduces a description of how some particular data element is stored and how many of them there are. Hence, in a file where there are 12 vertices, each represented as a floating point (X,Y,Z) triple, one would expect to see:
        element vertex 12
        property float x
        property float y
        property float z
        Other 'property' lines might indicate that colours or other data items are stored at each vertex and indicate the data type of that information. Regarding the data type there are two variants, depending on the source of the ply file. The type can be specified with one of char uchar short ushort int uint float double, or one of int8 uint8 int16 uint16 int32 uint32 float32 float64. For an object with ten polygonal faces, one might see:
        element face 10
        property list uchar int vertex_index
        The word 'list' indicates that the data is a list of values, the first of which is the number of entries in the list (represented as a 'uchar' in this case). In this example each list entry is represented as an 'int'. At the end of the header, there must always be the line:
        end_header
        https://en.wikipedia.org/wiki/PLY_(file_format)
    */

    // this method will return model =  {"v":vertices, "f":faces, "z":{"min":minz, "max":maxz}};
    static loadModel(content){
        var vertices = [];
        var faces = [];
        const allLines =  content.split(/\r\n|\n/);
        var vertexCount, faceCount;
        
        // Reading line by line
        var i = 0;
        if(allLines[i]!="ply"){
            console.log("Wrong data format, expect ply");
        }
        i+=1;
        while(allLines[i]!="end_header"){
            if(allLines[i].startsWith("format")){
                if(!allLines[i].startsWith("format ascii")){
                    console.log("We support ply ascii only")
                    return;
                }
            }
            if(allLines[i].startsWith("element vertex")){
                vertexCount = parseInt(allLines[i].split(" ")[2], 10)
            }
            if(allLines[i].startsWith("element face")){
                faceCount = parseInt(allLines[i].split(" ")[2], 10)
            }
            i+=1;
        }
        var vertexLineStart = i+1;
        
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
}