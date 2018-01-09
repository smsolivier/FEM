res = .1; 
res2 = .15; 
Point(1) = {0, 0, 0, res};
//+
Point(2) = {2, 0, 0, res};
//+
Point(3) = {2, 1, 0, res};
//+
Point(4) = {0, 1, 0, res};
//+
Point(5) = {1, 0.5, 0, 1.0};
//+
Point(6) = {.85, .5, 0, res};
//+
Point(7) = {1.15, .5, 0, res};
//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Circle(5) = {6, 5, 7};
//+
Circle(6) = {7, 5, 6};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Line Loop(2) = {6, 5};
//+
Plane Surface(1) = {1, 2};
//+
Physical Line(2) = {1};
//+
Physical Line(1) = {4, 2};
//+
Physical Line(3) = {3};
//+
Physical Line(1) += {6, 5};
//+
low = .125; 
high = .875; 
Physical Surface(7) = {1};
