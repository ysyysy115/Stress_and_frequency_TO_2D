function [dfreqfx]=Frequency_Sensitivity(ELE_COOR,ELE_MAP,KeStandard,MeStandard,x,nelx,nely,EIGVAL,TRUE_EIGVEC,p,q,weight,No_mode);

nel=nelx*nely;
dfreqfx=zeros(nely,nelx);
for mode=1:No_mode
% Assemble u
for iel=1:nel
    u(1,1,iel)=TRUE_EIGVEC(ELE_COOR(iel,1)*2-1,mode);
    u(2,1,iel)=TRUE_EIGVEC(ELE_COOR(iel,1)*2,mode);
    u(3,1,iel)=TRUE_EIGVEC(ELE_COOR(iel,2)*2-1,mode);
    u(4,1,iel)=TRUE_EIGVEC(ELE_COOR(iel,2)*2,mode);
    u(5,1,iel)=TRUE_EIGVEC(ELE_COOR(iel,3)*2-1,mode);
    u(6,1,iel)=TRUE_EIGVEC(ELE_COOR(iel,3)*2,mode);
    u(7,1,iel)=TRUE_EIGVEC(ELE_COOR(iel,4)*2-1,mode);
    u(8,1,iel)=TRUE_EIGVEC(ELE_COOR(iel,4)*2,mode);
end

% Sensitivity
for iel=1:nel
    [row,col]=find(ELE_MAP==iel);
    if x(row,col)>0.1
        dfreqfx1(iel)=u(:,:,iel)'*(p*x(row,col)^(p-1)*KeStandard-EIGVAL(mode,mode)*MeStandard)*u(:,:,iel);
    else
        dfreqfx1(iel)=u(:,:,iel)'*(p*x(row,col)^(p-1)*KeStandard-EIGVAL(mode,mode)*q*x(row,col)^(q-1)*MeStandard)*u(:,:,iel);
    end
end
dfreqfx1=dfreqfx1/max(dfreqfx1);
% dfreqfx1=dfreqfx1/1e7;

dfreqfx2=zeros(nely,nelx);
iel=0;
for col=1:nelx
    for row=1:nely
        iel=iel+1;
        dfreqfx2(row,col)=dfreqfx1(iel);
    end
end
dfreqfx2=flipud(dfreqfx2);

dfreqfx=weight(mode)*dfreqfx2+dfreqfx;

end
