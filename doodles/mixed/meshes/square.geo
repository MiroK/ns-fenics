Point(1) = {0., 0., 0., 1.};
Point(2) = {1., 0., 0., 1.};
Point(3) = {0., 1., 0., 1.};
Point(4) = {1., 1., 0., 1.};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};

Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Physical Line(7) = {4, 2}; // pressure
Physical Line(8) = {3, 1}; // wall, no slip
Physical Surface(9) = {6};
