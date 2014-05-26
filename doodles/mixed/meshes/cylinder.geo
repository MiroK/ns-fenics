// Gmsh project created on Wed May 21 14:44:51 2014

// Walls
l_wall = 1;
Point(1) = {0.0, 0.0, 0, l_wall};
Point(2) = {2.2, 0.0, 0, l_wall};
Point(3) = {2.2, 0.41, 0, l_wall};
Point(4) = {0.0, 0.41, 0, l_wall};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Cylinder
l_cyl = 1;
Point(5) = {0.2, 0.2, -0, l_cyl};
Point(6) = {0.15, 0.2, -0, l_cyl};
Point(7) = {0.25, 0.2, -0, l_cyl};
Point(8) = {0.2, 0.25, -0, l_cyl};
Point(9) = {0.2, 0.15, -0, l_cyl};


// Surface

Circle(5) = {6, 5, 8};
Circle(6) = {8, 5, 7};
Circle(7) = {7, 5, 9};
Circle(8) = {9, 5, 6};
Line Loop(9) = {4, 1, 2, 3};
Line Loop(10) = {7, 8, 5, 6};
Plane Surface(11) = {9, 10};
Physical Line(12) = {5, 6, 8, 7, 1, 3};
Physical Line(13) = {4};
Physical Surface(14) = {11};
Physical Line(15) = {2};
