% File: Q4_COARSE_MESH_DATA2
%
global nnd nel nne nelx nely nndx nndy nodof eldof n nf ngp lx ly NODE_MAP ELE_MAP ELE_COOR thick density
global geom connec dee nf Nodal_loads
%
% To change the size of the mesh, alter the next statements
%
nelx=60;
nely=20;
nndx=nelx+1;
nndy=nely+1;
lx=1; %length x of element
ly=1; %length y of element

nnd = nndx*nndy ; % Number of nodes:
nel = nelx*nely ; % Number of elements:
nne = 4 ; % Number of nodes per element:
nodof =2; % Number of degrees of freedom per node
ngp = 2 ;% number of Gauss points
eldof = nne*nodof; % Number of degrees of freedom per element

%
%% Generate Node Map
NODE_MAP=zeros(nndy,nndx);
inode=0;
for col=1:nndx
    for row=1:nndy
        inode=inode+1;
        NODE_MAP(row,col)=inode;
    end
end
NODE_MAP=flipud(NODE_MAP);

%% Generate Element Map and Element-node Coordinate(connectivity martix)
%node order
%43
%12
ELE_MAP=zeros(nely,nelx);
iel=0;
for col=1:nelx
    for row=1:nely
        iel=iel+1;
        ELE_MAP(row,col)=iel;
    end
end
ELE_MAP=flipud(ELE_MAP);

ELE_COOR=zeros(nel,4);
for iel=1:nel
    [hang,lie]=find(ELE_MAP==iel);
    ELE_COOR(iel,1)=NODE_MAP(hang+1,lie);
    ELE_COOR(iel,2)=NODE_MAP(hang+1,lie+1);
    ELE_COOR(iel,3)=NODE_MAP(hang,lie+1);
    ELE_COOR(iel,4)=NODE_MAP(hang,lie);
end

%% Nodes coordinates x and y
% geom = [0, -10.0; ... % x and y coordinates of node 1
% 0.0 0.0; ... % x and y coordinates of node 2
% 0.0 10.0; ... % x and y coordinates of node 3
% 10.0 -10.0; ... % x and y coordinates of node 4
% 10.0 0.0; ... % x and y coordinates of node 5
% 10.0 10.0; ... % x and y coordinates of node 6
% 20.0 -10.0; ... % x and y coordinates of node 7
% 20.0 0.0; ... % x and y coordinates of node 8
% 20.0 10.0; ... % x and y coordinates of node 9
% 30.0 -10.0; ... % x and y coordinates of node 10
% 30.0 0.0; ... % x and y coordinates of node 11
% 30.0 10.0; ... % x and y coordinates of node 12
% 40.0 -10.0; ... % x and y coordinates of node 13
% 40.0 0.0; ... % x and y coordinates of node 14
% 40.0 10.0; ... % x and y coordinates of node 15
% 50.0 -10.0; ... % x and y coordinates of node 16
% 50.0 0.0; ... % x and y coordinates of node 17
% 50.0 10.0; ... % x and y coordinates of node 18
% 60.0 -10.0; ... % x and y coordinates of node 19
% 60.0 0.0; ... % x and y coordinates of node 20
% 60.0 10.0]; % x and y coordinates of node 21

geom=zeros(nnd,2);
for inode=1:nnd
    [hang,lie]=find(NODE_MAP==inode);
    geom(inode,1)=(lie-1)*lx; %x coordinate
    geom(inode,2)=(nndy-hang)*ly; %y coordinate
end

%
%
%
% disp ('Nodes X-Y coordinates')
% geom
%
%% Element connectivity
% connec= [ 1 4 5 2 ;... % Element 1
% 2 5 6 3 ;... % Element 2
% 4 7 8 5 ;... % Element 3
% 5 8 9 6 ;... % Element 4
% 7 10 11 8 ;... % Element 5
% 8 11 12 9 ;... % Element 6
% 10 13 14 11 ;... % Element 7
% 11 14 15 12 ;... % Element 8
% 13 16 17 14 ;... % Element 9
% 14 17 18 15 ;... % Element 10
% 16 19 20 17 ;... % Element 11
% 17 20 21 18]; % Element 12

connec=ELE_COOR;

%
%
% disp ('Elements connectivity')
% connec
%
%% Material properties
E = 200000.; % Elastic modulus in MPa
vu = 0.3; % Poisson’s ratio
thick = 1.; % Beam thickness in mm
density = 8e-9 ; % Material Density in kg/m^3 or ton/mm^3
%
%% Form the elastic matrix for plane stress
%
dee = formdsig(E,vu);
%
%
%% Boundary conditions
%
nf = ones(nnd, nodof); % Initialize the matrix nf to 1
%Find the elements on the two lateral sides
left=NODE_MAP(:,1);
right=NODE_MAP(:,end);
for ind=1:size(left,1)
    nf(left(ind),:)=[0 0];
end
for ind=1:size(right,1)
    nf(right(ind),:)=[0 1];
end

% nf(17,:) = [0 0];
% nf(18,:) = [0 0];
% nf(19,:) = [0 0]; % Node 19 is restrained in the x and y directions
% nf(20,:) = [0 0]; % Node 20 is restrained in the x and y directions

%
%% Counting of the free degrees of freedom
%
n=0;
for i=1:nnd
for j=1:nodof
if nf(i,j) ~= 0
n=n+1;
nf(i,j)=n;
end
end
end
%
%% loading
%
Nodal_loads= zeros(nnd, 2); % Initialize the matrix of nodal loads to 0
%
% Apply a concentrated at the node having x = 0, and y = 0.
%
Force = 1000.; % N
%
Nodal_loads(1,:) = [0. -Force];
