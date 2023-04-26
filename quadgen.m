function [a,b,nodeTable,elementTable] = quadgen(Lx,Ly,ex,ey)

    % % % Input:
    % Lx = Length in x-direction
    % Ly = Length in y-direction
    % ex = Number of elements in x-direction
    % ey = Number of elements in y-direction
    % % % Output:
    % a = Element length in x-direction
    % b = Element length in y-direction
    % nodeTable    = Nodal coordinate table
    % elementTable = Element connectivity table
    % % %

    a = (Lx/ex)/2; b = (Ly/ey)/2;
    nx = ex+1; ny = ey+1;
    etot = ex*ey; ntot = nx*ny;
    
    temp = (1:ny:(nx-1)*ny)';
    temp = [temp temp+ny temp+ny+1 temp+1];
    temp = repmat(temp,ey,1) + kron([0:ey-1]',[ones(ex,4)]);
    elementTable = [(1:etot)' temp];
    
    xtemp = 0:2*a:Lx;
    xtemp = repmat(xtemp,ny,1);
    ytemp = (0:2*b:Ly)';
    ytemp = repmat(ytemp,1,nx);
    nodeTable = [(1:ntot)' xtemp(:) ytemp(:)];    
    
end