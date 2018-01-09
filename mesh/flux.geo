//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, .5, .5, 0};
//+
Line Loop(2) = {3, 4, 1, 2};
//+
Plane Surface(2) = {2};
//+
Physical Line(2) = {4, 1};
//+
Physical Line(1) = {3, 2};
//+
Physical Surface(5) = {1};
