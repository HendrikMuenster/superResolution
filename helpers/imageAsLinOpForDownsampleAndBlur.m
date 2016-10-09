%% A = imageAsLinOpForDownsampleAndBlur(kernelSize, img,fac)
% returns a matrix A such that downsample(k*img) by factor "fac"
% can be written as A vec(k), where k is a convolution kernel.
% This might help when optimizing for the kernel k
% to use for this procedure. k is of size (2*kernelSize+1)^2 if fac is odd
% and of size (2*kernelSize)^2 if fac is even to center the kernel
% symmetrically around the high resolution pixels that will result in the
% low resolution pixel. 

function A = imageAsLinOpForDownsampleAndBlur(kernelSize, img,fac)

if mod(fac,2)==1
    modifyShift = 0;
else
    modifyShift=-1;
end
additionalShift = ceil((fac-1)/2);
img = padarray(img, [kernelSize,kernelSize], 'replicate');
totalSize = 2*kernelSize+1+modifyShift;
%[ny,nx,nc]=size(img);
%A = zeros((1+nx-totalSize)*(1+ny-totalSize)*nc,totalSize);
for i=-kernelSize:(kernelSize+modifyShift)
    for j=-kernelSize:(kernelSize+modifyShift)
        temp = img((additionalShift+kernelSize+1+i):fac:(end+additionalShift-kernelSize+i+modifyShift), (additionalShift+kernelSize+1+j):fac:(end+additionalShift-kernelSize+j+modifyShift),:);
        ih = i + kernelSize+1;
        jt = j + kernelSize+1;
        A(:,ih+(jt-1)*totalSize) = temp(:);
    end
end
 
    
