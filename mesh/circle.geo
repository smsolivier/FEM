//+
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 1, 0, 2*Pi};
//+
Line Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Physical Line(1) = {1};
//+
Physical Surface(2) = {1};
