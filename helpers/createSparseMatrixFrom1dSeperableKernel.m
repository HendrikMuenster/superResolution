function K = createSparseMatrixFrom1dSeperableKernel(xDirKernel, yDirKernel, imDim1,imDim2, boundaryCondition)

if size(xDirKernel,2)==1
    xDirKernel=xDirKernel';
end

if size(yDirKernel,2)==1
    yDirKernel=yDirKernel';
end

nx = ceil((size(xDirKernel,2)-1)/2);
ny = ceil((size(yDirKernel,2)-1)/2);


%convolution of u with separable kernels is the same as Ay*u*Ax;
if strcmp(boundaryCondition,'Neumann')
    Ax = spdiags(repmat(xDirKernel, [imDim2,1]), 0:size(xDirKernel,2)-1, imDim2, imDim2+size(xDirKernel,2)-1)';
    Ax = [sum(Ax(1:nx+1,:),1);Ax(nx+2:end,:)];
    Ax = [Ax(1:end-nx-1,:);sum(Ax(end-nx:end,:),1)];
    Ay = spdiags(repmat(yDirKernel, [imDim1,1]), 0:size(yDirKernel,2)-1, imDim1, imDim1+size(yDirKernel,2)-1);
    Ay = [sum(Ay(:,1:ny+1),2),Ay(:,ny+2:end)];
    Ay = [Ay(:,1:end-ny-1),sum(Ay(:,end-ny:end),2)];
else
    Ax = spdiags(repmat(xDirKernel, [imDim2,2]), -nx:(size(xDirKernel,2)-1-nx), imDim2, imDim2)';
    Ay = spdiags(repmat(yDirKernel, [imDim1,2]), -ny:(size(yDirKernel,2)-1-ny), imDim1, imDim1);
end

K = kron(Ax',Ay);

