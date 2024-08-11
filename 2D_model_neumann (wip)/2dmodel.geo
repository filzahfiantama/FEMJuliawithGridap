// Coordinates for rectangle corners
Point(1) = {0, -3600, 0, 1.0};
Point(2) = {0, 0, 0, 1.0};
Point(3) = {7000, 0, 0, 1.0};
Point(4) = {7000, -3600, 0, 1.0};

// Lines connecting two points given on the right side
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Connecting lines, creating a full rectangle
Curve Loop(1) = {2, 3, 4, 1};

// Defining the plane, later labeled as domain
Plane Surface(1) = {1};
Physical Surface("domain", 5) = {1};

// Defining the boundaries and labeling it
Physical Curve("top", 6) = {2};
Physical Curve("sidewalls", 7) = {1, 3};
Physical Curve("bottom", 8) = {4};

// Preparing the geometry for meshing
/// Defining the density of the mesh
//// --> Dividing each side into n points, not how many rectangles!!!
Transfinite Curve {1, 2, 3, 4} = 10 Using Progression 1;
Transfinite Surface {1};
