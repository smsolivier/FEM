//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Line Loop(2) = {3, 4, 1, 2};
//+
Plane Surface(2) = {2};
//+
Physical Line(1) = {4, 3, 2, 1};
//+
Physical Surface(2) = {1};
