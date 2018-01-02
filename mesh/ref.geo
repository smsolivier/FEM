//+
ref1 = 1; 
ref2 = 1; 
ref3 = 1; 
Point(1) = {0, 0, 0, ref3};
//+
Point(2) = {1, 0, 0, ref3};
//+
Point(3) = {1, 1, 0, ref3};
//+
Point(4) = {0, 1, 0, ref3};
//+
Point(5) = {.3, .3, 0, ref2};
//+
Point(6) = {.7, .3, 0, ref2};
//+
Point(7) = {.7, .7, 0, ref2};
//+
Point(8) = {.3, .7, 0, ref2};
//+
Line(1) = {4, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 4};
//+
Line(5) = {8, 7};
//+
Line(6) = {7, 6};
//+
Line(7) = {6, 5};
//+
Line(8) = {5, 8};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Line Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(1) = {1, 2};
//+
Plane Surface(2) = {2};
//+
Physical Line("dirichlet") = {4, 1, 2, 3};
//+
//+
Physical Surface("core") = {1};
Physical Surface("reflector") = {2};

