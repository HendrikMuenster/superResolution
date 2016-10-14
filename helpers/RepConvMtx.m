function [K] = RepConvMtx(k,dimsLarge)
%hmm
%only square kernels with uneven size :<
Nx = dimsLarge(1); Ny = dimsLarge(2);

[kx,ky] = size(k);
klength = kx*ky;
K = spalloc(Nx*Ny,Nx*Ny,klength*Nx*Ny);
kex = floor(kx/2); %extent
kc  = ceil(kx/2);  % center


mpos = (1:Nx*Ny);
for i = -kex:kex %rows
    for j = -kex:kex %columns
        npos = max(min(mpos+i+j*Nx,Nx*Ny),1);
        tmpK  = sparse(mpos,npos,ones(Nx*Ny,1)*k(kc+i,kc+j),Nx*Ny,Nx*Ny);
        K     = K +tmpK;
    end
end


end

