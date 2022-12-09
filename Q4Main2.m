% THIS PROGRAM USES AN 4-NODDED QUADRILATERAL ELEMENT FOR THE LINEAR ELASTIC
% STATIC ANALYSIS OF A TWO DIMENSIONAL PROBLEM
%
% Make these variables global so they can be shared by other functions
%
function [Freq,TRUE_EIGVEC,EIGVAL,KeStandard,MeStandard,mm,kk]=Q4Main2(x,p,q,weight,No_mode)
global nnd nel nne nodof eldof n ngp
global geom connec dee nf Nodal_loads
%
format long g
%
% To change the size of the problem or change the elastic properties
% supply another input file
%
Q4_COARSE_MESH_DATA2
%
%%%%%%%%%%%%%%%%%%%%%%%%%% End of input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Assemble the global force vector
%
fg=zeros(n,1);
for i=1: nnd

if nf(i,1) ~= 0
fg(nf(i,1))= Nodal_loads(i,1);
end
if nf(i,2) ~= 0
fg(nf(i,2))= Nodal_loads(i,2);
end
end
%
% Form the matrix containing the abscissas and the weights of Gauss points
%
ngp = 2;
samp=gauss(ngp);
%
% Numerical integration and assembly of the global stiffness matrix
%
% initialize the global stiffness matrix to zero
kk = zeros(n, n);
mm = zeros(n, n);
%
for i=1:nel
[coord,g] = elem_q4(i) ; % coordinates of the nodes of element i,
% and its steering vector
ke=zeros(eldof,eldof) ; % Initialize the element stiffness matrix
% to zero
for ig=1: ngp
wi = samp(ig,2);
for jg=1: ngp
wj=samp(jg,2);
[der,fun] = fmlin(samp, ig,jg); % Derivative of shape functions
%in local coordinates
jac=der*coord; % Compute Jacobian matrix
d=det(jac); % Compute determinant of Jacobian
% matrix
jac1=inv(jac); % Compute inverse of the Jacobian
deriv=jac1*der; % Derivative of shape functions
% in global coordinates
bee=formbee(deriv,nne,eldof); % Form matrix [B]
ke=ke + d*thick*wi*wj*bee'*dee*bee; % Integrate stiffness matrix
end
end
KeStandard=ke;
me=q4_mass(lx,ly,thick,density);
% mdiag=repelem(2e-9,8);
% me=diag(mdiag);
MeStandard=me;

%
%
%SIMP
[row,col]=find(ELE_MAP==i);
ke=x(row,col)^p*ke;
if x(row,col)>0.1
    me=x(row,col)*me;
else
    me=x(row,col)^q*me;
end
k_local(:,:,iel)=ke;
m_local(:,:,iel)=me;
%
%
%

kk=form_kk(kk,ke, g); % assemble global stiffness matrix
mm=form_kk(mm,me, g); % assemble global mass matrix
end
%
%
%%%%%%%%%%%%%%%%%%%%%%% End of assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
delta = kk\fg ; % solve for unknown displacements
%
disp('node x_disp y_disp ') %
for i=1: nnd %
if nf(i,1) == 0 %
x_disp =0.; %
else
x_disp = delta(nf(i,1)); %
end
%
if nf(i,2) == 0 %
y_disp = 0.; %
else
y_disp = delta(nf(i,2)); %
end
% disp([i x_disp y_disp]) % Display displacements of each node
DISP(i,:) = [x_disp y_disp];
end
%
%
ngp=1; % Calculate stresses and strains at
%the center of each element
samp=gauss(ngp);
%
for i=1:nel
[coord,g] = elem_q4(i); % coordinates of the nodes of element i,
% and its steering vector
eld=zeros(eldof,1); % Initialize element displacement to zero
for m=1:eldof %
if g(m)==0 %
eld(m)=0.; %
else %
eld(m)=delta(g(m)); % Retrieve element displacement from the
% global displacement vector
end
end
%
for ig=1: ngp
wi = samp(ig,2);
for jg=1: ngp
wj=samp(jg,2);
[der,fun] = fmlin(samp, ig,jg); % Derivative of shape functions
% in local coordinates
jac=der*coord; % Compute Jacobian matrix
jac1=inv(jac); % Compute inverse of the Jacobian
deriv=jac1*der; % Derivative of shape functions
% in global coordinates
bee=formbee(deriv,nne,eldof); % Form matrix [B]
eps=bee*eld; % Compute strains
sigma=dee*eps; % Compute stresses
end
end
SIGMA(i,:)=sigma ; % Store stresses for all elements
end
%
% Average stresses at nodes
%
[ZX, ZY, ZT, Z1, Z2]=stresses_at_nodes_Q4(SIGMA);
%
%
% Plot stresses in the x_direction
%

% U2 = DISP(:,2);
% cmin = min(U2);
% cmax = max(U2);
% caxis([cmin cmax]);
% patch('Faces', connec, 'Vertices', geom, 'FaceVertexCData',U2,...
% 'Facecolor','interp','Marker','.');
% colorbar;

%% Vibration test
[EIGVEC, EIGVAL] = eigs(kk,mm,No_mode,'smallestabs'); %LANCZOS ALGORITHM FOR SMALLEST NINE EIGENVALS.

[DVEC, IND]      = sort(diag(EIGVAL));
EIGVAL           = EIGVAL(IND,IND);
EIGVEC           = EIGVEC(:,IND);
EIGVAL           = (1/(2*pi))*sqrt(EIGVAL);

% Mass Normalized EIGVEC
for mode=1:No_mode
    NF(mode)=sqrt(EIGVEC(:,mode)'*mm*EIGVEC(:,mode));
    EIGVEC(:,mode)=EIGVEC(:,mode)/NF(mode);
%     NF_CHECK(mode)=EIGVEC(:,mode)'*mm*EIGVEC(:,mode);
end

% Make up EIGVEC for BCs
[EIGVEC_node,EIGVEC_dof]=find(nf);
EIGLIST_TEMP=[EIGVEC_node,EIGVEC_dof];
EIGLIST_TEMP=sortrows(EIGLIST_TEMP,1);
EIGLIST=zeros(size(EIGVEC,1),1);
for i=1:size(EIGVEC,1)
    if EIGLIST_TEMP(i,2)==1;
        EIGLIST(i)=EIGLIST_TEMP(i,1)*2-1;
    else
        EIGLIST(i)=EIGLIST_TEMP(i,1)*2;
    end
end

TRUE_EIGVEC=zeros(2*nnd,No_mode);
for idof=1:2*nnd
    if ismember(idof,EIGLIST)
        TRUE_EIGVEC(idof,:)=EIGVEC(find(EIGLIST==idof),:);
    else
        TRUE_EIGVEC(idof,:)=zeros(1,No_mode);
    end
end

Freq=0;
for mode=1:No_mode
    Freq=Freq+weight(mode)*EIGVAL(mode,mode);
end

%%
%Plot
Scaled_EIGVEC=zeros(2*nnd,1);
geom_new=zeros(nnd,2);

figure(1)
for kk = 1:No_mode
subplot(3,3,kk)    

Scaled_EIGVEC=TRUE_EIGVEC(:,kk)/max(TRUE_EIGVEC(:,kk))*5;
%New Coor
for inode=1:nnd
    geom_new(inode,1)=geom(inode,1)+Scaled_EIGVEC(inode*2-1,1);
    geom_new(inode,2)=geom(inode,2)+Scaled_EIGVEC(inode*2,1);
%     level(inode,1)=Scaled_EIGVEC(inode*2-1,1)+Scaled_EIGVEC(inode*2,1);
end

scatter(geom_new(:,1),geom_new(:,2))

% shading interp
% colorbar
% xlabel('clamped edge / chord')
% ylabel('span length')
% zlabel('deflection w')
% view(45,45)
% axis equal
% daspect([1 1 3])
% rotate3d
% grid off
% title(['mode: ' num2str(kk) ' / '   num2str(EIGVAL(kk,kk)) ' Hz'],'FontSize',10)
% axis off
end

end


