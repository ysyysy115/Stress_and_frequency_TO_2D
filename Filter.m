function [H]=Filter(nel,nelx,nely,rmin,ELE_MAP,lx,ly)
    CENTROID_MAP_x=zeros(nel,1);
    CENTROID_MAP_y=zeros(nel,1);
    for iel=1:nel
        [hang,lie]=find(ELE_MAP==iel);
        gaodu=nely-hang+1;
        CENTROID_MAP_x(iel)=(lie-1)*lx+0.5*lx;
        CENTROID_MAP_y(iel)=(gaodu-1)*ly+0.5*ly;
    end
    H=zeros(nel,nel);
    for i=1:nel
        for j=1:nel
            distance=sqrt((CENTROID_MAP_x(i)-CENTROID_MAP_x(j))^2+(CENTROID_MAP_y(i)-CENTROID_MAP_y(j))^2);
            H(i,j)=max(0,rmin-distance);
        end
    end
