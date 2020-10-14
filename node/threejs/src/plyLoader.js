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
        var edge_set = [];
        for(var i=0; i<faceCount; i++){
            var line = allLines[i+faceLineStart];
            var infos = line.split(" ");
            v1 = parseInt(infos[1], 10);
            v2 = parseInt(infos[2], 10);
            v3 = parseInt(infos[3], 10);        
            faces.push([v1,v2,v3]);    
            // edges count
            var e1 = [v1,v2];
            var e2 = [v2,v3];
            var e3 = [v3,v1];
            var el = [e1,e2,e3];
            // check duplicate
            if(edge_set.length==0){
                edge_set.push(e1);
            }                        
            for(var j=0;j<el.length;j++){
                var exists = false;
                for(var k=0;k<edge_set.length;k++){
                    if(el[j][0]==edge_set[k][0]&&el[j][1]==edge_set[k][1]){
                        exists = true;
                    }
                    if(el[j][0]==edge_set[k][1]&&el[j][1]==edge_set[k][0]){
                        exists = true;
                    }        
                }
                if(!exists){
                    edge_set.push(el[j])
                }
            }
        }
        var model =  {"v":vertices, "f":faces, "e":edge_set, "z":{"min":minz, "max":maxz}};
        console.log(model);
        return model;
    }

    static exportModel(model){
        var res = "ply\n";
        res += "comment Model exported by WebTerrain\n";
        res += "format ascii 1.0\n";
        res += "element vertex " + vertices.length + "\n";
        res += "element face " + faces.length + "\n";
        res += "end_header\n";
        var vertices = model.v;
        var faces = model.f;
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