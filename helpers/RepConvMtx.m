function [K] = RepConvMtx(k,dimsLarge)
% hmm
% only kernels with uneven size :<
%
%
if sum(mod(size(k),2)) ~=2
    error('kernel has even dimensions (which is not implemented)')
end

Nx = dimsLarge(1); Ny = dimsLarge(2);
N = Nx*Ny;

[kx,ky] = size(k);
klength = kx*ky;
kex = floor(kx/2); key = floor(ky/2); % extent
kcx  = ceil(kx/2); kcy = ceil(ky/2);  % center
%evenx = 1-mod(kx,2); %??
%eveny = 1-mod(ky,2);

SparseAlloc = zeros(N*klength,3);
aC = 1;              % allocation counter
mpos = (1:N)';
south_border   = ceil(mpos/Nx)*Nx;  %>= mpos
north_border   = floor((mpos-1)/Nx)*Nx+1;  %<= mpos


for i = -kex:kex     % rows
    for j = -key:key % columns
        npos  = max(min(mpos+i,south_border),north_border); %all of these are vectors
        west_border    = mod(npos-1,Nx)+1; 
        east_border    = mod(npos-1,Nx)+1+Nx*(Ny-1);
        npos = max(min(npos+j*Nx,east_border),west_border);
        SparseAlloc(N*(aC-1)+1:N*aC,:) = [mpos,npos,ones(N,1)*k(kcx-i,kcy-j)];
        aC = aC +1;
    end
end

K = sparse(SparseAlloc(:,1),SparseAlloc(:,2),SparseAlloc(:,3),N,N); %Sum and place into sparse matrix


end

