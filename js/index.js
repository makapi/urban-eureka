var container, width, height, n;
var mesh, camera, scene, renderer, controls, trackball, ambientLight, directionalLight;
var colors = [];
var u, v, dx, dt, count;
var h, s, l;
var initialPosition;

window.addEventListener("load", function(){
    init();
    animate();
});

function init(){
    container = document.createElement('div');
    document.body.appendChild(container);
    width = window.innerWidth;
    height = window.innerHeight;
    n = 128;

    init_scene();
    init_object();
    init_camera();
    init_light();
    init_renderer();

    function init_scene(){
        scene = new THREE.Scene();
    }

    function init_camera(){
        camera = new THREE.PerspectiveCamera(45, width / height, 1, 10000);
        camera.position.set(1, 0, 6);
        scene.add(camera);
        trackball = new THREE.TrackballControls(camera);
        trackball.noRotate = false; 
        trackball.rotateSpeed = 0.3;
        trackball.noZoom = false; 
        trackball.zoomSpeed = 0.3;
        trackball.noPan = false; 
        trackball.panSpeed = 0.3;
        trackball.staticMoving = true; 
        trackball.dynamicDampingFactor = 0.3;
    }

    function init_light(){
        directionalLight = new THREE.DirectionalLight(0xaaaaaa, 0.6);
        directionalLight.position.set(1, 1, 1);
        scene.add(directionalLight);
        directionalLight2 = new THREE.DirectionalLight(0xaaaaaa, 0.6);
        directionalLight2.position.set(-1, -1, -1);
        scene.add(directionalLight2);
        ambientLight = new THREE.AmbientLight(0x999999);
        scene.add(ambientLight);
    }

    function init_renderer(){
        renderer = new THREE.WebGLRenderer({antialias:true, alpha:true});
        renderer.setClearColor(0xeeeeee);
        renderer.setSize(width, height);
        container.appendChild(renderer.domElement);
    }

    function init_object(){
        mesh = mesh2();
        scene.add(mesh);
        
        //var axis = new THREE.AxisHelper(1000);
        //scene.add(axis);

        //var cShere = new CustomSphere();
        //scene.add(cShere);
    }

    window.addEventListener('resize', onWindowResize, false);
}

function animate(){
    requestAnimationFrame(animate);
    render();
    updateMesh2();
    trackball.update();
}

function render(){
    renderer.clear();
    renderer.render(scene, camera);   
}

function onWindowResize(){
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();  
    renderer.setSize(window.innerWidth, window.innerHeight);
}

function mesh2(){
    var geometry, material, _mesh;

    dt = 0.001;
    dx = 0.01;
    u = new Array(3);
    v = [];
    count = 0;

    h = 0.5;
    s = 1;
    l = 1;

    initialPosition = [];

    function mapN(k, imin, imax, omin, omax){
        return ((k - imin) * (omax - omin) / (imax - imin) + omin);
    }

    geometry = new THREE.Geometry();

    u[0] = [];
    u[1] = [];
    u[2] = [];

    for(var i=0; i<n; i++){
        u[0][i] = [];
        u[1][i] = [];
        u[2][i] = [];
        v[i] = new Array(n);
        initialPosition[i] = [];

        for(var j=0; j<n; j++){ 
            var ry = mapN(i, 0, n-1, Math.PI/2, 3 * Math.PI /2);
            var rx = mapN(j, 0, n-1, 0, 2 * Math.PI);
            var x = Math.cos(ry) * Math.sin(rx);
            var z = Math.cos(ry) * Math.cos(rx);
            var y = Math.sin(ry);

            var pos = new THREE.Vector3(x, y, z);
            geometry.vertices.push(pos);
            colors.push(new THREE.Color().setHSL(h, s, l));

            u[0][i][j] = new THREE.Vector3(pos.x, pos.y, pos.z);
            u[1][i][j] = new THREE.Vector3(pos.x, pos.y, pos.z);
            u[2][i][j] = new THREE.Vector3(pos.x, pos.y, pos.z);
            v[i][j] = 1;   
            initialPosition[i][j] = new THREE.Vector3();
            initialPosition[i][j].x = pos.x;
            initialPosition[i][j].y = pos.y;
            initialPosition[i][j].z = pos.z;
        }
    }
    
    for(var i=0; i<n-1; i++){
        for(var j=0; j<n-1; j++){ 
            var color1 = [];
            var color2 = [];
            color1[0] = colors[(n*i)+j];
            color1[1] = colors[(n*i)+(j+1)];
            color1[2] = colors[(n*(i+1))+j];
            color2[0] = colors[(n*i)+(j+1)];
            color2[1] = colors[(n*(i+1))+(j+1)];
            color2[2] = colors[(n*(i+1))+j];
            geometry.faces.push(new THREE.Face3((n*i)+j, (n*i)+(j+1), (n*(i+1))+j, null, color1));
            geometry.faces.push(new THREE.Face3((n*i)+(j+1), (n*(i+1))+(j+1), (n*(i+1))+j, null, color2));
        }
    }
    for(var i=0; i<n-2; i++){
        var j = n-1;
        var color1 = [];
        var color2 = [];
        color1[0] = colors[(n*i)+j];
        color1[1] = colors[(n*i)];
        color1[2] = colors[(n*(i+1))+j];
        color2[0] = colors[(n*i)];
        color2[1] = colors[(n*(i+1))];
        color2[2] = colors[(n*(i+1))+j];
        geometry.faces.push(new THREE.Face3((n*i)+j, (n*i), (n*(i+1))+j, null, color1));
        geometry.faces.push(new THREE.Face3((n*i), (n*(i+1)), (n*(i+1))+j, null, color2));
    }

    geometry.computeFaceNormals();
    geometry.computeVertexNormals();

    material = new THREE.MeshPhongMaterial({//color: 0xafeeee, 
                                            vertexColors: THREE.VertexColors,
                                            ambient: 0xafeeee, 
                                            //ambient: 0x000000,
                                            side: THREE.DoubleSide, 
                                            //specular: 0x080808, 
                                            specular: 0xffffff, 
                                            shininess: 250});
    _mesh = new THREE.Mesh(geometry, material);
    return _mesh;
}

var nz1 = Math.random()*Math.random();
var nz2 = Math.random()*Math.random();
var nz3 = Math.random()*Math.random();
var nz4 = Math.random()*Math.random();
var nz5 = Math.random()*Math.random();
var nz6 = Math.random()*Math.random();


function updateMesh2(){      
    function gaussian(x, y, z, mx, my, mz){
        //var sigma = 0.01;
        var sigma = 0.05;
        var g = Math.exp(-1*((x-mx)*(x-mx)+(y-my)*(y-my)+(z-mz)*(z-mz))/(2*sigma))/Math.sqrt(2*Math.PI*sigma);
        return g*0.004;
    }

    function unitVec(x, y, z){
        var mag = Math.sqrt(x*x + y*y + z*z);
        var u_v = new THREE.Vector3(x/mag, y/mag, z/mag);
        return u_v;
    }

    for(var k=0; k<10; k++){
        count++;
        if((count%8)==0){
            var offset = 16;

            nz1 += (0.01);
            nz2 += (0.02);
            nz3 += (0.03);
            nz4 += 0.01;
            nz5 += 0.01;
            nz6 += 0.01;
            var mi = Math.round(((noise.perlin3(nz1,nz2,nz3)+1)/2)*(n-1));
            var mj = Math.round(((noise.perlin3(nz2,nz3,nz1)+1)/2)*(n-1));
            //mi %= (n-1);
            //mj %= (n-1);

            if(mi<offset){
                for(var i=(mi-offset); i<(mi+offset); i++){
                    for(var j=(mj-offset); j<(mj+offset); j++){
                        var ii, ij;
                        if(i<0){
                            ii = i*(-1);
                            ij = j+(n/2);
                            if(ij>(n-1))
                                ij%=n;
                            var gs = gaussian(u[1][ii][ij].x, u[1][ii][ij].y, u[1][ii][ij].z, u[1][mi][mj].x, u[1][mi][mj].y, u[1][mi][mj].z); 
                            u[2][ii][ij].x += initialPosition[ii][ij].x *gs;
                            u[2][ii][ij].y += initialPosition[ii][ij].y *gs;
                            u[2][ii][ij].z += initialPosition[ii][ij].z *gs;
                        }  
                        else{
                            ii = i;
                            if(j<0)
                                ij = n+j;
                            else if(j>(n-1))
                                ij = j%n;
                            else
                                ij = j;
                            var gs = gaussian(u[1][ii][ij].x, u[1][ii][ij].y, u[1][ii][ij].z, u[1][mi][mj].x, u[1][mi][mj].y, u[1][mi][mj].z);
                            
                            u[2][ii][ij].x += initialPosition[ii][ij].x *gs;
                            u[2][ii][ij].y += initialPosition[ii][ij].y *gs;
                            u[2][ii][ij].z += initialPosition[ii][ij].z *gs;
                        }
                    }
                }
            }
            else if(mi>(n-1)-offset){
                for(var i=(mi-offset); i<(mi+offset); i++){
                    for(var j=(mj-offset); j<(mj+offset); j++){
                        var ij, jj;
                        if(i>(n-1)){
                            ii = i-(n-1);
                            ij = j+(n/2);
                            if(ij>(n-1))
                                ij%=n;
                            var gs = gaussian(u[1][ii][ij].x, u[1][ii][ij].y, u[1][ii][ij].z, u[1][mi][mj].x, u[1][mi][mj].y, u[1][mi][mj].z);
                            
                            u[2][ii][ij].x += initialPosition[ii][ij].x *gs;
                            u[2][ii][ij].y += initialPosition[ii][ij].y *gs;
                            u[2][ii][ij].z += initialPosition[ii][ij].z *gs;
                        }
                        else{
                            ii = i;
                            if(j<0)
                                ij = n+j;
                            else if(j>(n-1))
                                ij = j%n;
                            else
                                ij = j;
                            var gs = gaussian(u[1][ii][ij].x, u[1][ii][ij].y, u[1][ii][ij].z, u[1][mi][mj].x, u[1][mi][mj].y, u[1][mi][mj].z);
                            
                            u[2][ii][ij].x += initialPosition[ii][ij].x *gs;
                            u[2][ii][ij].y += initialPosition[ii][ij].y *gs;
                            u[2][ii][ij].z += initialPosition[ii][ij].z *gs;

                        }
                    }
                }
            }
            else{
                for(var i=(mi-offset); i<(mi+offset); i++){
                    for(var j=(mj-offset); j<(mj+offset); j++){
                        var ij;
                        if(j<0)
                            ij = n + j;
                        else if(j>(n-1))
                            ij = j % n;
                        else
                            ij = j;

                        var gs = gaussian(u[1][i][ij].x, u[1][i][ij].y, u[1][i][ij].z, u[1][mi][mj].x, u[1][mi][mj].y, u[1][mi][mj].z);
                       
                        u[2][i][ij].x += initialPosition[i][ij].x *gs;
                        u[2][i][ij].y += initialPosition[i][ij].y *gs;
                        u[2][i][ij].z += initialPosition[i][ij].z *gs;
                    }
                }
            }
        }
        else{
            for(var i=0; i<n; i++){
                for(var j=0; j<n; j++){
                    if(i==0){
                        if(j==0){
                            u[2][i][j].x = 2.0 * u[1][i][j].x - u[0][i][j].x 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i+1][j].x + u[1][i][j+1].x + u[1][i][n-1].x - 4.0 * u[1][i][j].x);

                            u[2][i][j].y = 2.0 * u[1][i][j].y - u[0][i][j].y 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i+1][j].y + (1) + u[1][i][j+1].y + u[1][i][n-1].y - 4.0 * u[1][i][j].y);

                            u[2][i][j].z = 2.0 * u[1][i][j].z - u[0][i][j].z 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i+1][j].z + u[1][i][j+1].z + u[1][i][n-1].z - 4.0 * u[1][i][j].z);

                        }
                        else if(j==n-1){

                            u[2][i][j].x = 2.0 * u[1][i][j].x - u[0][i][j].x 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i+1][j].x + u[1][i][0].x + u[1][i][j-1].x - 4.0 * u[1][i][j].x);

                            u[2][i][j].y = 2.0 * u[1][i][j].y - u[0][i][j].y 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i+1][j].y +(1) + u[1][i][0].y + u[1][i][j-1].y - 4.0 * u[1][i][j].y);

                            u[2][i][j].z = 2.0 * u[1][i][j].z - u[0][i][j].z 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i+1][j].z + u[1][i][0].z + u[1][i][j-1].z - 4.0 * u[1][i][j].z);

                        }
                        else{

                            u[2][i][j].x = 2.0 * u[1][i][j].x - u[0][i][j].x 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i+1][j].x + u[1][i][j+1].x + u[1][i][j-1].x - 4.0 * u[1][i][j].x);

                            u[2][i][j].y = 2.0 * u[1][i][j].y - u[0][i][j].y 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i+1][j].y + (1) + u[1][i][j+1].y + u[1][i][j-1].y - 4.0 * u[1][i][j].y);

                            u[2][i][j].z = 2.0 * u[1][i][j].z - u[0][i][j].z 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i+1][j].z + u[1][i][j+1].z + u[1][i][j-1].z - 4.0 * u[1][i][j].z);

                        }

                    }
                    else if(i==n-1){   
                        if(j==0){

                            u[2][i][j].x = 2.0 * u[1][i][j].x - u[0][i][j].x 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i-1][j].x + u[1][i][j+1].x + u[1][i][n-1].x - 4.0 * u[1][i][j].x);

                            u[2][i][j].y = 2.0 * u[1][i][j].y - u[0][i][j].y 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i-1][j].y + (-1) + u[1][i][j+1].y + u[1][i][n-1].y - 4.0 * u[1][i][j].y);

                            u[2][i][j].z = 2.0 * u[1][i][j].z - u[0][i][j].z 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i-1][j].z + u[1][i][j+1].z + u[1][i][n-1].z - 4.0 * u[1][i][j].z);

                        }
                        else if(j==n-1){

                            u[2][i][j].x = 2.0 * u[1][i][j].x - u[0][i][j].x 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i-1][j].x + u[1][i][0].x + u[1][i][j-1].x - 4.0 * u[1][i][j].x);

                            u[2][i][j].y = 2.0 * u[1][i][j].y - u[0][i][j].y 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i-1][j].y + (-1) + u[1][i][0].y + u[1][i][j-1].y - 4.0 * u[1][i][j].y);

                            u[2][i][j].z = 2.0 * u[1][i][j].z - u[0][i][j].z 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i-1][j].z + u[1][i][0].z + u[1][i][j-1].z - 4.0 * u[1][i][j].z);

                        }
                        else{

                            u[2][i][j].x = 2.0 * u[1][i][j].x - u[0][i][j].x 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i-1][j].x + u[1][i][j+1].x + u[1][i][j-1].x - 4.0 * u[1][i][j].x);

                            u[2][i][j].y = 2.0 * u[1][i][j].y - u[0][i][j].y 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i-1][j].y + (-1) + u[1][i][j+1].y + u[1][i][j-1].y - 4.0 * u[1][i][j].y);

                            u[2][i][j].z = 2.0 * u[1][i][j].z - u[0][i][j].z 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i-1][j].z + u[1][i][j+1].z + u[1][i][j-1].z - 4.0 * u[1][i][j].z);

                        }
                    }
                    else if(j==0){
                        u[2][i][j].x = 2.0 * u[1][i][j].x - u[0][i][j].x 
                               + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                               * (u[1][i+1][j].x + u[1][i-1][j].x + u[1][i][j+1].x + u[1][i][n-1].x - 4.0 * u[1][i][j].x);

                        u[2][i][j].y = 2.0 * u[1][i][j].y - u[0][i][j].y 
                               + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                               * (u[1][i+1][j].y + u[1][i-1][j].y + u[1][i][j+1].y + u[1][i][n-1].y - 4.0 * u[1][i][j].y);

                        u[2][i][j].z = 2.0 * u[1][i][j].z - u[0][i][j].z 
                               + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                               * (u[1][i+1][j].z + u[1][i-1][j].z + u[1][i][j+1].z + u[1][i][n-1].z - 4.0 * u[1][i][j].z); 


                    }
                    else if(j==n-1){
                        u[2][i][j].x = 2.0 * u[1][i][j].x - u[0][i][j].x 
                               + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                               * (u[1][i+1][j].x + u[1][i-1][j].x + u[1][i][0].x + u[1][i][j-1].x - 4.0 * u[1][i][j].x);

                        u[2][i][j].y = 2.0 * u[1][i][j].y - u[0][i][j].y 
                               + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                               * (u[1][i+1][j].y + u[1][i-1][j].y + u[1][i][0].y + u[1][i][j-1].y - 4.0 * u[1][i][j].y);

                        u[2][i][j].z = 2.0 * u[1][i][j].z - u[0][i][j].z 
                               + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                               * (u[1][i+1][j].z + u[1][i-1][j].z + u[1][i][0].z + u[1][i][j-1].z - 4.0 * u[1][i][j].z);
                    }
                    else{
                        u[2][i][j].x = 2.0 * u[1][i][j].x - u[0][i][j].x 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i+1][j].x + u[1][i-1][j].x + u[1][i][j+1].x + u[1][i][j-1].x - 4.0 * u[1][i][j].x);

                        u[2][i][j].y = 2.0 * u[1][i][j].y - u[0][i][j].y 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i+1][j].y + u[1][i-1][j].y + u[1][i][j+1].y + u[1][i][j-1].y - 4.0 * u[1][i][j].y);

                        u[2][i][j].z = 2.0 * u[1][i][j].z - u[0][i][j].z 
                                   + v[i][j] * v[i][j] * dt * dt / (dx * dx) 
                                   * (u[1][i+1][j].z + u[1][i-1][j].z + u[1][i][j+1].z + u[1][i][j-1].z - 4.0 * u[1][i][j].z);
                    }
                }
            }
        }
        
        if(count%300==0){
            for(var i=0; i<n; i++){
                for(var j=0; j<n; j++){
                    u[2][i][j].x += (initialPosition[i][j].x - u[2][i][j].x)/90;
                    u[2][i][j].y += (initialPosition[i][j].y - u[2][i][j].y)/90;
                    u[2][i][j].z += (initialPosition[i][j].z - u[2][i][j].z)/90;  
                }
            }
        }
        
        for (var i=0; i<n; i++) {
            for (var j=0; j<n; j++){
                u[0][i][j].x = u[1][i][j].x;
                u[0][i][j].y = u[1][i][j].y;
                u[0][i][j].z = u[1][i][j].z;
                u[1][i][j].x = u[2][i][j].x;
                u[1][i][j].y = u[2][i][j].y;
                u[1][i][j].z = u[2][i][j].z;
            }
        }

    }
    

    var a = 0;
    for(var i=0; i<n; i++){
        for(var j=0; j<n; j++){
                var x = u[1][i][j].x;
                var y = u[1][i][j].y;
                var z = u[1][i][j].z;

                mesh.geometry.vertices[a].x = x;
                mesh.geometry.vertices[a].y = y;
                mesh.geometry.vertices[a].z = z;

                a++;           
        }
    }
    var f = 0;
    for(var i=0; i<n-1; i++){
        for(var j=0; j<n-1; j++){
            var x = mesh.geometry.vertices[(i*n)+j].x;
            var y = mesh.geometry.vertices[(i*n)+j].y;
            var z = mesh.geometry.vertices[(i*n)+j].z;
            var c = Math.sqrt(x*x+y*y+z*z);
            h = c / 3.5;
            //l = 0.8 + (Math.sqrt(c)-1);
            //l = 1-(Math.sqrt(c)-1);
            l = 0.91;
            mesh.geometry.faces[f].vertexColors[0] = new THREE.Color().setHSL(h,s,l);
            mesh.geometry.faces[f].vertexColors[1] = new THREE.Color().setHSL(h,s,l);
            mesh.geometry.faces[f].vertexColors[2] = new THREE.Color().setHSL(h,s,l);
            f++;
            mesh.geometry.faces[f].vertexColors[0] = new THREE.Color().setHSL(h,s,l);
            mesh.geometry.faces[f].vertexColors[1] = new THREE.Color().setHSL(h,s,l);
            mesh.geometry.faces[f].vertexColors[2] = new THREE.Color().setHSL(h,s,l);
            f++;
        }
    }

    for(var i=0; i<n-2; i++){
        var j = n-1;
        var x = mesh.geometry.vertices[(i*n)+j].x;
        var y = mesh.geometry.vertices[(i*n)+j].y;
        var z = mesh.geometry.vertices[(i*n)+j].z;
        var c = Math.sqrt(x*x+y*y+z*z);
        h = c / 3.5;
        //l = 0.8 + (Math.sqrt(c)-1);
        //l = 1-(Math.sqrt(c)-1);
        l = 0.91;
        mesh.geometry.faces[f].vertexColors[0] = new THREE.Color().setHSL(h,s,l);
        mesh.geometry.faces[f].vertexColors[1] = new THREE.Color().setHSL(h,s,l);
        mesh.geometry.faces[f].vertexColors[2] = new THREE.Color().setHSL(h,s,l);
        f++;
        mesh.geometry.faces[f].vertexColors[0] = new THREE.Color().setHSL(h,s,l);
        mesh.geometry.faces[f].vertexColors[1] = new THREE.Color().setHSL(h,s,l);
        mesh.geometry.faces[f].vertexColors[2] = new THREE.Color().setHSL(h,s,l);
        f++;
    }

    mesh.geometry.colorsNeedUpdate = true;

    mesh.geometry.normalsNeedUpdate = true;
    mesh.geometry.verticesNeedUpdate = true;

    mesh.geometry.computeFaceNormals();
    mesh.geometry.computeVertexNormals();   
}
