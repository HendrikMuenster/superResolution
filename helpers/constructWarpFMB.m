function [warpingOp] = constructWarpFMB(v)
%
% Create forward vs minus backward warp as seen in the old jointSuperResolution class
%



% Get sizes
[Ny,Nx,~,n4] = size(v);
nc =  n4+1;

spX = []; spY = []; spAlloc = [];
for i=1:nc - 1
    %extract flow field between i and i+1
    singleField = squeeze(v(:,:,i,:));
    
    %create warping operator forward and backward
    warpF = warpingOperator([Nx,Ny],0.5*singleField);
    warpB = warpingOperator([Nx,Ny],-0.5*singleField);
    
    %find out of range warps in each of the operators and set the
    %corresponding line in the other operator also to zero
    idF = sum(abs(warpF),2) == 0;
    idB = sum(abs(warpB),2) == 0;
    
    marker = idF + idB;
    marker = marker > 0;
    
    warpF(marker,:) = 0;
    warpB(marker,:) = 0;
    
    
    % Find nonzeros and add
    [sp_x,sp_y,walloc] = find(warpB);
    spX = [spX;(Nx*Ny*(i-1)+1)+sp_x-1]; %#ok<*AGROW>
    spY = [spY;(Nx*Ny*(i-1)+1)+sp_y-1];
    spAlloc = [spAlloc;walloc];
    
    [sp_x,sp_y,walloc] = find(warpF);
    spX = [spX;(Nx*Ny*(i-1)+1)+sp_x-1];
    spY = [spY;(Nx*Ny*i+1)+sp_y-1];
    spAlloc = [spAlloc;-walloc];
    
end
% Construct operator
warpingOp = sparse(spX,spY,spAlloc,Nx*Ny*nc,Nx*Ny*nc);
end

