function  spU = imageStackOp(u,kernelsize)
% IMAGESTACKOP
% Computes the appropriate matrix U, so that Downsample(U*vec(k)) = f
% u is a stack of grayscale images, [Nx,Ny,nc]
% k is a vector of vectorized blur kernels k_i

% size(U) = [Nx*Ny*nc,kernelsize^2*nc]

% written with uneven kernel sizes in mind 
[Ny,Nx,nc] = size(u);
padsize = floor(kernelsize/2);
klength = kernelsize^2;

StackOp = spalloc(Nx*Ny*nc,klength*nc,Nx*Ny*klength*nc);
for i = 1:nc
    uP     = padarray(u(:,:,i),[padsize,padsize],'replicate');
    tempOp = im2col(uP,[kernelsize,kernelsize],'sliding')';
    StackOp((Nx*Ny*(i-1)+1):(Nx*Ny*i),(klength*(i-1)+1):(klength*i)) = sparse(tempOp); %#ok<SPRIX>
end
spU = StackOp;
end