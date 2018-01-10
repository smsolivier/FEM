//+
SetFactory("OpenCASCADE");
Circle(1) = {0, -0, 0, 1, 0, 2*Pi};
//+
Circle(2) = {0, -0, 0, 1.1, 0, 2*Pi};
//+
Circle(3) = {0, 0, 0, 1.25, 0, 2*Pi};
//+
Rectangle(1) = {-2, -2, 0, 4, 4, 0};
//+
Line Loop(2) = {1};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {2};
//+
Line Loop(4) = {1};
//+
Plane Surface(3) = {3, 4};
//+
Line Loop(5) = {3};
//+
Line Loop(6) = {2};
//+
Plane Surface(4) = {5, 6};
//+
Line Loop(7) = {6, 7, 4, 5};
//+
Line Loop(8) = {3};
//+
Plane Surface(5) = {7, 8};
//+
Physical Line("dirichlet",2) = {4,5,6,7}; 
//Physical Line("dirichlet") = {4, 6};
//Physical Line("reflecting") = {5,7}; 
Physical Surface("core") = {2}; 
Physical Surface("gap") = {3}; 
Physical Surface("clad") = {4}; 
Physical Surface("moderator") = {5}; 
Characteristic Length { 0 } = .1;
