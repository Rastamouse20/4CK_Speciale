function [filter] = filterPrep(nelx,nely,nodes,rmin)
    iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
    jH = ones(size(iH));
    sH = zeros(size(iH));
        
    k = 0;
    for i1 = 1:nelx
      for j1 = 1:nely
        e1 = (j1-1)*nelx+i1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
          for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
            e2 = (j2-1)*nelx+i2;
            k = k+1;
            iH(k) = e1;
            jH(k) = e2;
            sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
          end
        end
      end
    end

    filter.H = sparse(iH,jH,sH);
    filter.Hs = sum(filter.H,2);
end