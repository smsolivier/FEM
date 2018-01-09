res = .075; 
Point(1) = {0, 0, 0, .1};
//+
Point(2) = {.5, 0, 0, .1};
//+
Point(3) = {.5, .4, 0, res};
//+
Point(4) = {.5, .5, 0, res};
//+
Point(5) = {.4, .5, 0, res};
//+
Point(6) = {.1, .5, 0, res};
//+
Point(7) = {0, .5, 0, res};
//+
Point(8) = {0, .4, 0, res};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 1};
//+
Line Loop(1) = {5, 6, 7, 8, 1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Line(1) = {7, 8, 1, 2, 3};
//+
Physical Line(2) = {5, 6, 4};
//+
Physical Surface(5) = {1};
