// this code was taken from J&J's github; generates ellipsoidal mesh for FEM cell sim.

SetFactory("OpenCASCADE");

l = 100; // Side length of box
Sphere(1) = {0, 0, 0, 1};

// principal semi-major axes
psax = 11.5;
psay = 7.6;
psaz = 18.75;
Dilate {{0, 0, 0}, {psax, psay, psaz}} {Volume{1};}

Sphere(100) = {0, 0, 0, 1};

Box(2) = {-l/2, -l/2, -l/2, l, l, l};

//BooleanDifference(30) = {Volume{2}; Delete; }{Volume{1}; Delete; }; //deletes the cytoplasm; I want to keep it.
//BooleanDifference(3) = {Volume{1}; Delete; }{Volume{100}; Delete; }; //tries to delete nucleus; not 100% sure if needed.

// I don't get why this is working...physical volumes are not showing up.

Physical Volume(100) = {2}; //the gel; with everything.
Physical Volume(300) = {100}; //nucleus; I think I want it to be empty though - why are there errors when empty, and not with
//john and jp's cell??
Physical Volume(200) = {1}; //ellipsoidal cell (cytoplasm)

Physical Surface(201) = {7};
Physical Surface(200) = {4, 5, 3, 2, 6, 1};

Mesh.Algorithm = 6;
Characteristic Length{:} = 10;
// giving different characteristic length factors to different surfaces/volumes (cool!)
Characteristic Length{PointsOf{Physical Volume{4};}}  = 10;
Characteristic Length{PointsOf{Surface{1};}} = 1; //meshing the box layer finely which is weird, was 7 i changed to 1.

// Generate Mesh
//Mesh 3;
//Mesh.MshFileVersion = 2.2;
//Save "ellipsoid.msh";




