function [me]=q4_mass(lx,ly,thickness,density)
%This function computes the element mass matrix of a linear quadralateral
%element

me=density*lx*ly*thickness/36*...
   [4 0 2 0 1 0 2 0;...
    0 4 0 2 0 1 0 2;...
    2 0 4 0 2 0 1 0;...
    0 2 0 4 0 2 0 1;...
    1 0 2 0 4 0 2 0;...
    0 1 0 2 0 4 0 2;...
    2 0 1 0 2 0 4 0;...
    0 2 0 1 0 2 0 4];

end