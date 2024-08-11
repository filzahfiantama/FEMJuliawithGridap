// Coordinates for rectangle corners
Point(1) = {0, -0, 0, 1.0};
Point(2) = {-0, 1, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {1, -0, 0, 1.0};

// Lines connecting two points given on the right side
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Connecting lines, creating a full rectangle
Curve Loop(1) = {2, 3, 4, 1};

// Defining the plane, later labeled as domain
Plane Surface(1) = {1};
Physical Surface("domain", 6) = {1};

// Defining the boundary and labeling it
Physical Curve("boundary", 5) = {1, 2, 3, 4};

// Preparing the geometry for meshing
Transfinite Surface {1} = {2, 3, 4, 1};
// Defining the density of the mesh
// --> Dividing each side into n points, not how many rectangles!!!
Transfinite Curve {1, 3, 2, 4} = 9 Using Progression 1;
