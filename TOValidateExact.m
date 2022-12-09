%% Dimensions
lx=1;
ly=1;
No_el_x=500;
No_el_y=45;
No_el=No_el_x*No_el_y;
No_nd_x=No_el_x+1;
No_nd_y=No_el_y+1;
No_nd=No_nd_x*No_nd_y;
threshold=0.75;

%% Preprocessing
ELE_MAP=1:No_el;
ELE_MAP=transpose(reshape(ELE_MAP,[No_el_x,No_el_y]));
NODE_MAP=1:No_nd;
NODE_MAP=transpose(reshape(NODE_MAP,[No_nd_x,No_nd_y]));

ELE_COOR=zeros(No_el_y,No_el_x,4);
ELE_COOR(:,:,1)=NODE_MAP(1:No_el_y,1:No_el_x);
ELE_COOR(:,:,2)=ELE_COOR(:,:,1)+1;
ELE_COOR(:,:,3)=NODE_MAP(2:No_el_y+1,1:No_el_x);
ELE_COOR(:,:,4)=ELE_COOR(:,:,3)+1;

active_el=[];
SolidFlag=zeros(No_el_y,No_el_x);
for iel=1:No_el
    [row,col]=find(ELE_MAP==iel);
    if x(row,col)>=threshold
        SolidFlag(row,col)=1;
        active_el=[active_el,iel];
    end
end

left_el=[];
right_el=[];
left_nd=[];
right_nd=[];
for ilayer=1:No_el_y
    if SolidFlag(ilayer,1)
        left_el=[left_el,ELE_MAP(ilayer,1)];
    end
    if SolidFlag(ilayer,end)
        right_el=[right_el,ELE_MAP(ilayer,end)];
    end
end
for iel=1:size(left_el,2)
    [row,col]=find(ELE_MAP==left_el(iel));
    left_nd=[left_nd,ELE_COOR(row,col,1),ELE_COOR(row,col,3)];
end
for iel=1:size(right_el,2)
    [row,col]=find(ELE_MAP==right_el(iel));
    right_nd=[right_nd,ELE_COOR(row,col,2),ELE_COOR(row,col,4)];
end
left_nd=unique(left_nd);
right_nd=unique(right_nd);

active_nd=[];
for iel=1:size(active_el,2)
    [row,col]=find(ELE_MAP==active_el(iel));
    active_nd=[active_nd,ELE_COOR(row,col,1),ELE_COOR(row,col,2),ELE_COOR(row,col,3),ELE_COOR(row,col,4)];
end
active_nd=unique(active_nd);

%% Write input file
fid = fopen('D:\Matlabtemp\TOvalidateExact.inp','w');

fprintf(fid,'*HEADING\n');
fprintf(fid,'*Preprint, echo=NO, model=NO, history=NO, contact=NO\n');

%% Generate Nodes
fprintf(fid,'*Node\n');
for ilayer=1:No_nd_y
    fprintf(fid,'%i,0,%f,0\n',NODE_MAP(ly,1),-(ilayer-1)*ly);
    fprintf(fid,'%i,%f,%f,0\n',NODE_MAP(ly,end),No_el_x*lx,-(ilayer-1)*ly);
end
for ilayer=1:No_nd_y
    fprintf(fid,'*NGEN\n');
    fprintf(fid,'%i,%i,1\n',NODE_MAP(ly,1),NODE_MAP(ly,end));
end

%% Generate Elements
No_Solid=0;
for iel=1:No_el
    [row,col]=find(ELE_MAP==iel);
    if x(row,col)>=threshold
        No_Solid=No_Solid+1;
        fprintf(fid,'*Element, type=S4, ELSET=SolidEL\n');
        fprintf(fid,'%i,%i,%i,%i,%i\n',iel,ELE_COOR(row,col,1),ELE_COOR(row,col,2),ELE_COOR(row,col,4),ELE_COOR(row,col,3));
    end
end

%% Sets
fprintf(fid,'*Nset, nset=left\n');
for inode=1:size(left_nd,2)
    fprintf(fid,'%i,\n',left_nd(inode));
end

fprintf(fid,'*Nset, nset=right\n');
for inode=1:size(right_nd,2)
    fprintf(fid,'%i,\n',right_nd(inode));
end

fprintf(fid,'*Nset, nset=lateral\n');
for inode=1:size(active_nd,2)
    fprintf(fid,'%i,\n',active_nd(inode));
end

%% Generate sections
fprintf(fid,'*Shell Section, elset=SolidEL, material=Material-1\n');
fprintf(fid,'%f, 5\n',5); %thickness, # of integration points

fprintf(fid,'*Material, name=Material-1\n');
fprintf(fid,'*Density\n');
fprintf(fid,'7.85e-9,\n');
fprintf(fid,'*Elastic\n');
fprintf(fid,'210000, 0.3\n');

%% Step 2 Frequency
fprintf(fid,'*Step, name=Step-2, nlgeom=NO, perturbation\n');
fprintf(fid,'*Frequency, eigensolver=Subspace, sim, acoustic coupling=off, normalization=mass\n');
fprintf(fid,'4, , , , , \n'); % # of modes required

% BCs
fprintf(fid,'*Boundary\n');
% Clamped left
fprintf(fid,'left, 1, 6\n');
% No rotations of lateral sides
fprintf(fid,'lateral, 3, 3\n');
% fprintf(fid,'lateral, 4, 4\n');
% fprintf(fid,'lateral, 5, 5\n');
% fprintf(fid,'lateral, 6, 6\n');
% No x, no rx right
fprintf(fid,'right, 1, 1\n');
fprintf(fid,'right, 6, 6\n');

fprintf(fid,'*Restart, write, frequency=0\n');
fprintf(fid,'*Output, field, variable=PRESELECT\n');
fprintf(fid,'*End Step\n');

fclose(fid);
disp('!!! SUCCESSFULLY WRITING INPUT FILE !!!');

!del D:\temp\TOvalidateExact.*
!copy D:\Matlabtemp\TOvalidateExact.inp D:\temp\TOvalidateExact.inp
cd('D:\temp')
!G:\SIMULIA\Commands\abaqus job=TOvalidateExact
