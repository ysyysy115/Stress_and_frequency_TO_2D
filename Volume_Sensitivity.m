function [fval,dfdx]=Volume_Sensitivity(volfrac,x,nelx,nely,H,Hs)
dv = ones(nely,nelx)/(nelx*nely);
dv(:) = H*(dv(:)./Hs);
fval=[mean(x(:))-volfrac];
dfdx=[dv(:)'];
