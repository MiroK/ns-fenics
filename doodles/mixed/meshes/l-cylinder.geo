// Gmsh project created on Wed May 21 14:44:51 2014

u = 0.1;
pi = 3.141592653589793;
s = 2;
t = s - 4*u*Sin(pi/4);

Point(1) = {4*u*Cos(5*pi/4), 4*u*Sin(5*pi/4), 0, 1}; // P
Point(2) = {4*u*Cos(6*pi/4), 4*u*Sin(6*pi/4), 0, 1}; // D
Point(3) = {4*u*Cos(7*pi/4), 4*u*Sin(7*pi/4), 0, 1}; // R
Point(4) = {4*u*Cos(2*pi/4), 4*u*Sin(2*pi/4), 0, 1}; // V

Point(5) = {4*u*Cos(5*pi/4) + s*Sin(5*pi/4), 4*u*Sin(5*pi/4) - s*Cos(5*pi/4), 0, 1}; // Q
Point(6) = {4*u*Cos(7*pi/4) - s*Sin(7*pi/4), 4*u*Sin(7*pi/4) + s*Cos(7*pi/4), 0, 1}; // S

Point(7) = {4*u*Cos(2*pi/4) + t*Sin(5*pi/4), 4*u*Sin(2*pi/4) - t*Cos(5*pi/4), 0, 1}; // Q
Point(8) = {4*u*Cos(2*pi/4) - t*Sin(7*pi/4), 4*u*Sin(2*pi/4) + t*Cos(7*pi/4), 0, 1}; // V
 

Point(9) = {0., 0., 0, 1.0};
Point(10) = {2*u, 0., 0, 1.0};
Point(11) = {-2*u, 0., 0, 1.0};
Point(12) = {0., 2*u, 0, 1.0};
Point(13) = {0., -2*u, 0, 1.0};

Circle(1) = {12, 9, 10};
Circle(2) = {10, 9, 13};
Circle(3) = {13, 9, 11};
Circle(4) = {11, 9, 12};
Line(5) = {5, 1};
Line(6) = {7, 4};
Line(7) = {4, 8};
Line(8) = {6, 3};
Line(9) = {8, 6};
Line(10) = {7, 5};

Circle(11) = {1, 9, 3};
Line Loop(12) = {6, 7, 9, 8, -11, -5, -10};
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(14) = {12, 13};
