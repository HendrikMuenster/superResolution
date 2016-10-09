function A = writeKernelToSparseDownsamplingMatrix(ker,fac,nx,ny)
%totKerSize = size(ker,1);
startX = 1+floor(fac/2);
startY = 1+floor(fac/2);

[x,y]=meshgrid(startX:fac:nx, startY:fac:ny);
nr =1;
if mod(fac,2)==1
    correction = 0;
else 
    correction=1;
end
temp = ceil((size(ker,1)-1)/2);

for i=-temp:(temp-correction)
    for j=-temp:(temp-correction)
        yCoord = min(max(y(:)+i,1),ny);
        xCoord = min(max(x(:)+j,1),nx);
        indi2(:,nr) = yCoord(:)+(xCoord(:)-1)*ny;
        nr = nr + 1;
    end
end
rows = size(indi2,1);
indi = repmat((1:rows)', [1,numel(ker)]);
vals = repmat(ker(:)', [rows,1]);

A = sparse(indi(:),indi2(:),vals(:),rows,nx*ny);
end
