L = .5;
res = .05; 
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {L, 0, 0, 1.0};
//+
Point(3) = {L, L, 0, res};
//+
Point(4) = {0, L, 0, res};
//+
Line(1) = {4, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 4};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Line(1) = {4, 3, 2};
//+
Physical Line(2) = {1};
//+
Physical Surface(5) = {1};
