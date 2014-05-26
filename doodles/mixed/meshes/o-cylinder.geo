u = 0.1;
pi = 3.141592653589793;

// Circle points
Point(1) = {0, 0, 0, 1};
Point(2) = {2.5*u*Cos(2*pi/4), 2.5*u*Sin(2*pi/4), 0, 1};
Point(3) = {2.5*u*Cos(4*pi/4), 2.5*u*Sin(4*pi/4), 0, 1};
Point(4) = {2.5*u*Cos(6*pi/4), 2.5*u*Sin(6*pi/4), 0, 1};
Point(5) = {2.5*u*Cos(8*pi/4), 2.5*u*Sin(8*pi/4), 0, 1};

// Lower turn points
Point(6) = {4*u*Cos(1*pi/4), 4*u*Sin(1*pi/4), 0, 1};
Point(7) = {4*u*Cos(2*pi/4), 4*u*Sin(2*pi/4), 0, 1};
Point(8) = {4*u*Cos(3*pi/4), 4*u*Sin(3*pi/4), 0, 1};

// Upper turn points
Point(9) = {4*u*Cos(5*pi/4), 4*u*Sin(5*pi/4), 0, 1};
Point(10) = {4*u*Cos(6*pi/4), 4*u*Sin(6*pi/4), 0, 1};
Point(11) = {4*u*Cos(7*pi/4), 4*u*Sin(7*pi/4), 0, 1};

// Lower left arm point
s = 3;
Translate {s*Sin(5*pi/4), -s*Cos(5*pi/4), 0 } { Duplicata{ Point{9}; } }

// Upper left arm point
t = s - 4*u;
Translate {t*Sin(5*pi/4), -t*Cos(5*pi/4), 0 } { Duplicata{ Point{8}; } }

// Lower right arm point
S = 3;
Translate {-S*Sin(7*pi/4), S*Cos(7*pi/4), 0 } { Duplicata{ Point{11}; } }

// Upper right arm point
T = S - 4*u;
Translate {-T*Sin(7*pi/4), T*Cos(7*pi/4), 0 } { Duplicata{ Point{6}; } }
Coherence;
Line(1) = {12, 9};
Line(2) = {13, 8};
Line(3) = {6, 15};
Line(4) = {14, 11};
Line(5) = {15, 14};
Line(6) = {13, 12};
Circle(7) = {9, 1, 11};
Circle(8) = {8, 1, 6};
Circle(9) = {3, 1, 2};
Circle(10) = {2, 1, 5};
Circle(11) = {5, 1, 4};
Circle(12) = {4, 1, 3};
Line Loop(13) = {2, 8, 3, 5, 4, -7, -1, -6};
Line Loop(14) = {9, 10, 11, 12};
Plane Surface(15) = {13, 14};
Physical Surface(16) = {15};
Physical Line(17) = {1, 2, 8, 3, 4, 7, 11, 10, 12, 9};
Physical Line(18) = {6};
Physical Line(19) = {5};
