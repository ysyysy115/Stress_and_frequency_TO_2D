% Frequency constrained topology optimization of LASEnelx=80*reso;
% Jason Yan, Jasmin Jelvica, 2022/06/23
% Linear quadratic elements are used

clc
clear all

global nnd nel nne nelx nely nndx nndy nodof eldof nv ngp lx ly NODE_MAP ELE_MAP ELE_COOR thick density No_mode
global geom connec dee nf Nodal_loads

reso=1;
nelx=50*reso;
nely=15*reso;
volfrac=0.3;
rmin=2;
p=3;
q=7;
No_mode=1;

nel=nelx*nely;
nele=nel;

Q4_COARSE_MESH_DATA2

%% Filter
H = Filter(nel,nelx,nely,rmin,ELE_MAP,lx,ly);
Hs = sum(H,2);

%% INITIALIZE
maxloop = 1000;    % Maximum number of iterations
tolx = 0.001;      % Terminarion criterion
displayflag = 0;  % Display structure flag

x(1:nely,1:nelx) = volfrac; 
change=1;
loop=0;

% INITIALIZE MMA OPTIMIZER
m     = 1;                % The number of general constraints.
nv     = nele;             % The number of design variables x_j.
xmin  = zeros(nv,1);       % Column vector with the lower bounds for the variables x_j.
xmax  = ones(nv,1);        % Column vector with the upper bounds for the variables x_j.
xold1 = x(:);             % xval, one iteration ago (provided that iter>1).
xold2 = x(:);             % xval, two iterations ago (provided that iter>2).
low   = ones(nv,1);        % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = ones(nv,1);        % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                % The constants a_0 in the term a_0*z.
a     = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 10000*ones(m,1);  % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.

while change > tolx  
  loop = loop + 1;
  xold = x;

 x(1,:)=1;
x(15,:)=1;

%WEIGHT SETTING
% weight(1:9)=1/9;
weight(1)=1;
weight(2)=1;
weight(3)=1;
weight(4)=1;
weight(5)=1;
weight(6)=1;
weight(7)=1;
weight(8)=1;
weight(9)=1;

xtemp=flipud(x);
xval  = xtemp(:);

% FE-ANALYSIS
  [Freq,TRUE_EIGVEC,EIGVAL,KeStandard,MeStandard,mm,kk]=Q4Main2(x,p,q,weight,No_mode);       
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
      c = -Freq;
      [dfreqfx]=Frequency_Sensitivity(ELE_COOR,ELE_MAP,KeStandard,MeStandard,x,nelx,nely,EIGVAL,TRUE_EIGVEC,p,q,weight,No_mode);
      dc = -dfreqfx;
% FILTERING OF SENSITIVITIES
  [dc]   = check(nelx,nely,rmin,x,dc);

f0val=c;
 [dc]   = check(nelx,nely,rmin,x,dc);
 dc=flipud(dc);
df0dx=dc(:);

% df0dx(:) = H*(df0dx(:)./Hs);

% VOLUME FRACTION CONSTRAINT AND SENSITIVITY
    [fval,dfdx]=Volume_Sensitivity(volfrac,x,nelx,nely,H,Hs);

    [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low,upp] = ...
mmasub(m, nv, loop, xval, xmin, xmax, xold1, xold2, ...
f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
% Update MMA Variables
xnew     = reshape(xmma,nely,nelx);
xold2    = xold1(:);
xold1    = x(:);

    change = max(abs(xnew(:)-x(:)));
    x = xnew;

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


