% Frequency constrained topology optimization of LASEnelx=80*reso;
% Jason Yan, Jasmin Jelvica, 2022/06/23
% Linear quadratic elements are used

clc
clear all

global nnd nel nne nelx nely nndx nndy nodof eldof n ngp lx ly NODE_MAP ELE_MAP ELE_COOR thick density No_mode
global geom connec dee nf Nodal_loads

reso=1;
nelx=60*reso;
nely=20*reso;
volfrac=0.3;
rmin=1.5;
p=3;
q=6;
No_mode=4;

nel=nelx*nely;

Q4_COARSE_MESH_DATA2

%% INITIALIZE
% INITIAL GUESS
% x(1:nely,1:nelx) = volfrac; %evenly distributed
% [x,v] = randfixedsum(nel,1,volfrac*nel,0,1); %randomly distributed
% x=reshape(x,nely,nelx);
x(1:nely,1:nelx) = volfrac; 
change=1;
loop=0;
maxloop=600;

while change > 0.01 && loop<maxloop 
  loop = loop + 1;
  xold = x;

x(1,:)=1;
x(20,:)=1;
%WEIGHT SETTING
% weight(1:9)=1/9;
weight(1)=1;
weight(2)=0.2;
weight(3)=0.2;
weight(4)=0.2;
weight(5)=1;
weight(6)=1;
weight(7)=1;
weight(8)=1;
weight(9)=1;

% FE-ANALYSIS
  [Freq,TRUE_EIGVEC,EIGVAL,KeStandard,MeStandard,mm,kk]=Q4Main2(x,p,q,weight,No_mode);       
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
      c = -Freq;
      [dfreqfx]=Frequency_Sensitivity(ELE_COOR,ELE_MAP,KeStandard,MeStandard,x,nelx,nely,EIGVAL,TRUE_EIGVEC,p,q,weight,No_mode);
      dc = -dfreqfx;
% FILTERING OF SENSITIVITIES
  [dc]   = check(nelx,nely,rmin,x,dc);    
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
  [x]    = OC(nelx,nely,x,volfrac,dc); 
% PRINT RESULTS
  change = max(max(abs(x-xold)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])
% PLOT DENSITIES  
figure(2)
  colormap(gray); imagesc(-real(x)); axis equal; axis tight; axis off;pause(1e-6);
end

% save data
% saveas(gcf,'result.fig');


