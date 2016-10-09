function [A, newKer] = optimalBlurAndDownsamplingKernel(obj,kernelSize, k)

A = imageAsLinOpForDownsampleAndBlur(kernelSize,obj.u(:,:,k),obj.factor);
temp = obj.imageSequenceSmall(:,:,k); 
newKer = reshape(lsqnonneg(A,temp(:)),[2*kernelSize+mod(kernelSize,2),2*kernelSize+mod(kernelSize,2)]);
%newKer = reshape((A'*A)\A'*temp(:), [2*kernelSize+mod(kernelSize,2),2*kernelSize+mod(kernelSize,2)]);
%newKer = max(newKer,0);
%newKer = newKer./sum(sum(newKer));
%sum((A*newKer(:) - temp(:)).^2)

[ny,nx] = size(obj.u(:,:,k));
A = writeKernelToSparseDownsamplingMatrix(newKer,obj.factor,nx,ny);
