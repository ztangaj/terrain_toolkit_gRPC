// after load model, we have vertices, can render the points
function onLoad() {
    console.log("onload");
    if(highlightVertex){				
        console.log(vertices);
        var colors = new Float32Array( vertices.length );
        var sizes = new Float32Array( vertices.length );

        var color = new THREE.Color();
        for ( var i = 0, l = vertices.length / 3; i < l; i ++ ) {
            // color.setHSL( 0.01 + 0.1 * ( i / l ), 1.0, 0.5 );
            color = new THREE.Color("rgb(255,0,0)");
            color.toArray( colors, i * 3 );
            sizes[ i ] = PARTICLE_SIZE;
        }
        var geometry = new THREE.BufferGeometry();
        geometry.addAttribute( 'position', new THREE.BufferAttribute( vertices, 3 ) );
        geometry.addAttribute( 'customColor', new THREE.BufferAttribute( colors, 3 ) );
        geometry.addAttribute( 'size', new THREE.BufferAttribute( sizes, 1 ) );

        // deprecated, don't use sprit for vertex now
        particles = new THREE.Points( geometry, pointMaterial );
        particles.name = "points";
        scene.add( particles );
        refreshControlPanel(panel_options);
    }
};

function onProgress( xhr ) {
    if ( xhr.lengthComputable ) {
        var percentComplete = xhr.loaded / xhr.total * 100;
        console.log( 'model ' + Math.round( percentComplete, 2 ) + '% downloaded' );
    }
}

function onError() {}		

var loader = new THREE.OBJLoader();					

loader.load( MODEL_PATH, function( object ) {	
    var boundingbox = new THREE.Box3();
    boundingbox.setFromObject(object);
    var max_z = boundingbox.max.z;
    var min_z = boundingbox.min.z;
    object.traverse( function ( child ) {
        if(child.isMesh) {
            child.material = generateMaterial(min_z, max_z);							
            vertices = child.geometry.attributes.position.array;
        };
    } );
    scene.add(object);									

    center = boundingbox.getCenter();
    distance = getCameraDistance(object);

    camera.position.set( center.x , center.y, center.z + distance);	
    controls.target.set(center.x, center.y, center.z );
    controls.update();						
}
, onLoad, onProgress, onError );				

