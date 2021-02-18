// got this method to make an ellipsoid in gmsh from forums.

//testing iteratively to figure out how to make ellipsoid with nucleus
//don't want the nucleus to have a mesh on the inside; need to verify this.

SetFactory("OpenCASCADE");

l = 50; // Side length of box
lcar3 = 2; // to comly with cheeshole format

Sphere(1) = {0, 0, 0, 10}; // was .5 ; I think I need to check out the  sizing of the mesh...not sure how to do this.
Dilate {{0, 0, 0}, {1.5, 1, 0.5}} { Volume{1}; }

Box(2) = {-l/2, -l/2, -l/2, l, l, l};
Sphere(3) = {0, 0, 0, 1}; //nucleus

// this line should remove the space that would be the nucleus
BooleanDifference(10) = { Volume{1}; Delete; }{ Volume{3}; Delete; }; //{ Volume{3}; Delete; }; //not sure if this will work
Surface Loop(11) = {10}; //not sure if this will work


//might want to include the gel as 100 too. also needs to conform to the code we used...oof.
//Physical Volume(300) = {3}; //nucleus
//Physical Volume(200) = {10}; //the ellipsoid with a hole in the center.
//likely want to have a nucleus with some boundary condition that the nucleus doesn't move, or is harder to move, idk. need to figure this out.


//not sure why there are 2 separate surfaces for the cube
//Physical Surface(201) = {7};
//Physical Surface(200) = {4, 5, 3, 2, 6, 1};

Volume(31) = {3};
Volume(21) = {11};

Physical Surface(201) = {7};
Physical Surface(101) = {4, 5, 3, 2, 6, 1}; //no clue how this works
Physical Surface(202) = {10}; //not sure if this will work - doesn't seem to work. should run my old code and see how it works.

Physical Volume(300) = {31};  //nucleus
Physical Volume(200) = {21}; //cell (ellipsoid with hole in center)
Physical Volume(100) = {1};  //gel

// mesh size determination
Characteristic Length{:} = 1; //0.1;
Characteristic Length{PointsOf{Physical Volume{1};}}  = 10;
//Characteristic Length{PointsOf{Surface{202};}} = 1;
Mesh.Algorithm = 6;

// this works/an ellipsoid is made. just need to make a cell surface too I think;
// not sure if the ellipsoid has a surface as of now.