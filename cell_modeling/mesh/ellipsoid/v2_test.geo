//need to figure out how to center the cell in the frame, and how to place the nucleus in the cell properly.

//Merge "cytod_uncentered_unpca.stl";
//CreateTopology;

//creates physical surface
//Physical Surface ( expression | char-expression <, expression> ) <+|->= { expression-list };
//"expression" = tag for surface
//"expression list" = tags of elementary surfaces that need to be grouped in the pysical surface

//Physical Surface(1) = {1}; //box representing gel domain

SetFactory("OpenCASCADE");

l = 50; // Side length of box
lcar3 = 2; // to comly with cheeshole format

Sphere(1000) = {0, 0, 0, 10}; // was .5 ; I think I need to check out the  sizing of the mesh...not sure how to do this.
Physical Surface(1) = {1000}; //IF I CAN MAKE THE ELLIPSOID PHYSICAL SURFACE 1, I SHOULD BE GOOD TO GO.
Dilate {{0, 0, 0}, {1.5, 1, 0.5}} { Volume{1000}; }

SetFactory("Built-in");

//not sure why it's not working with opencascade, might need to change things up aka generate ellipsoid differently.


//specify dimensions of the gel

length = 150.;
height = 150.;
depth = 150.;
lcar1 = 0.5; //was 8 then 0.5 - does it matter what this is??? didn't make a difference when i changed it.
lcar3 = 2; //was 32; I doubled each other them to make the mesh less fine.

//instantiate 8 verticies of the gel (cube)

//Point(newp) = {length,height,depth,lcar1}; /* Point      1 */
//Point(newp) = {length,height,0,lcar1}; /* Point      2 */
//Point(newp) = {0,height,depth,lcar1}; /* Point      3 */
//Point(newp) = {0,0,depth,lcar1}; /* Point      4 */
//Point(newp) = {length,0,depth,lcar1}; /* Point      5 */
//Point(newp) = {length,0,0,lcar1}; /* Point      6 */
//Point(newp) = {0,height,0,lcar1}; /* Point      7 */
//Point(newp) = {0,0,0,lcar1}; /* Point      8 */

//updated verticies of the gel/cube to be centered around the cell, which is located at the origin. don't know how to move it.
Point(newp) = {length/2,height/2,depth/2,lcar1}; /* Point      1 */
Point(newp) = {length/2,height/2,-depth/2,lcar1}; /* Point      2 */
Point(newp) = {-length/2,height/2,depth/2,lcar1}; /* Point      3 */
Point(newp) = {-length/2,-height/2,depth/2,lcar1}; /* Point      4 */
Point(newp) = {length/2,-height/2,depth/2,lcar1}; /* Point      5 */
Point(newp) = {length/2,-height/2,-depth/2,lcar1}; /* Point      6 */
Point(newp) = {-length/2,height/2,-depth/2,lcar1}; /* Point      7 */
Point(newp) = {-length/2,-height/2,-depth/2,lcar1}; /* Point      8 */

//connect these 8 points with 12 lines to form cube edges

Line(1) = {3,1};
Line(2) = {3,7};
Line(3) = {7,2};
Line(4) = {2,1};
Line(5) = {1,5};
Line(6) = {5,4};
Line(7) = {4,8};
Line(8) = {8,6};
Line(9) = {6,5};
Line(10) = {6,2};
Line(11) = {3,4};
Line(12) = {8,7};

//fill in each surface of the cube

Line Loop(13) = {-6,-5,-1,11};
Plane Surface(14) = {13};
Line Loop(15) = {4,5,-9,10};
Plane Surface(16) = {15};
Line Loop(17) = {-3,-12,8,10};
Plane Surface(18) = {17};
Line Loop(19) = {7,12,-2,11};
Plane Surface(20) = {19};
Line Loop(21) = {-4,-3,-2,1};
Plane Surface(22) = {21};
Line Loop(23) = {8,9,6,7};
Plane Surface(24) = {23};

Macro CheeseHole

  // In the following commands we use the reserved variable name `newp', which
  // automatically selects a new point tag. Analogously to `newp', the special
  // variables `newl', `newll, `news', `newsl' and `newv' select new curve,
  // curve loop, surface, surface loop and volume tags.
  //
  // If `Geometry.OldNewReg' is set to 0, the new tags are chosen as the highest
  // current tag for each category (points, curves, curve loops, ...), plus
  // one. By default, for backward compatibility, `Geometry.OldNewReg' is set
  // to 1, and only two categories are used: one for points and one for the
  // rest.

  p1 = newp; Point(p1) = {x,  y,  z,  lcar3};
  p2 = newp; Point(p2) = {x+r,y,  z,  lcar3};
  p3 = newp; Point(p3) = {x,  y+r,z,  lcar3};
  p4 = newp; Point(p4) = {x,  y,  z+r,lcar3};
  p5 = newp; Point(p5) = {x-r,y,  z,  lcar3};
  p6 = newp; Point(p6) = {x,  y-r,z,  lcar3};
  p7 = newp; Point(p7) = {x,  y,  z-r,lcar3};

  c1 = newc; Circle(c1) = {p2,p1,p7}; c2 = newc; Circle(c2) = {p7,p1,p5};
  c3 = newc; Circle(c3) = {p5,p1,p4}; c4 = newc; Circle(c4) = {p4,p1,p2};
  c5 = newc; Circle(c5) = {p2,p1,p3}; c6 = newc; Circle(c6) = {p3,p1,p5};
  c7 = newc; Circle(c7) = {p5,p1,p6}; c8 = newc; Circle(c8) = {p6,p1,p2};
  c9 = newc; Circle(c9) = {p7,p1,p3}; c10 = newc; Circle(c10) = {p3,p1,p4};
  c11 = newc; Circle(c11) = {p4,p1,p6}; c12 = newc; Circle(c12) = {p6,p1,p7};

  l1 = newll; Curve Loop(l1) = {c5,c10,c4};
  l2 = newll; Curve Loop(l2) = {c9,-c5,c1};
  l3 = newll; Curve Loop(l3) = {c12,-c8,-c1};
  l4 = newll; Curve Loop(l4) = {c8,-c4,c11};
  l5 = newll; Curve Loop(l5) = {-c10,c6,c3};
  l6 = newll; Curve Loop(l6) = {-c11,-c3,c7};
  l7 = newll; Curve Loop(l7) = {-c2,-c7,-c12};
  l8 = newll; Curve Loop(l8) = {-c6,-c9,c2};

  s1 = news; Surface(s1) = {l1};
  s2 = news; Surface(s2) = {l2};
  s3 = news; Surface(s3) = {l3};
  s4 = news; Surface(s4) = {l4};
  s5 = news; Surface(s5) = {l5};
  s6 = news; Surface(s6) = {l6};
  s7 = news; Surface(s7) = {l7};
  s8 = news; Surface(s8) = {l8};

  // We then store the surface loops tags in a list for later reference (we will
  // need these to define the final volume):

  theloops[t] = newsl;
  Surface Loop(theloops[t]) = {s1, s2, s3, s4, s5, s6, s7, s8};

  thehole = newv;
  Volume(thehole) = theloops[t];

Return

//x,y,z: (x,y,z) coordinate of the center of the nucleus
//r: the radius of the nucleus
//t: index

// ORIGINAL X,Y,Z,R,T VALUES
//x = 80; y = 73; z = 40; r = 3; t = 1;

x = 0; y = 0; z = 0; r = 3; t = 1;

Call CheeseHole;

Surface Loop(25) = {1,14,24,-18,22,16,-20};
Volume(26) = {25};
Surface Loop(27) = {-1,-s1,-s2,-s3,-s4,-s5,-s6,-s7,-s8}; 
Volume(28) = {27};

//not entirely sure where the specifc cell geometry comes from

Physical Surface(101) = {14,16,18,20,24};
Physical Surface(201) = {1};
Physical Volume(100) = {26}; // gel
Physical Volume(200) = {28}; // cytoplasm
Physical Volume(300) = thehole; // nucleus

Mesh.CharacteristicLengthFactor = 800; //100 = coarse, .5 was original (is very dense)...is this working is the question