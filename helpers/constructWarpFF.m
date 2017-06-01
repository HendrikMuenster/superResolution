function warpingOp = constructWarpFF(v,type,occ)
%Create warp operator for FF scheme
%type = 'F-I';
%type = 'I-F';
% v is of size Ny x Nx x numFrames-1 x 2 
% occ is of size Ny x Nx x numFrames-1

% Get sizes
[Ny,Nx,n4,~] = size(v);
nc =  n4+1;


if nargin < 3
    occ = [];
end


% create warping operator from given v
spX = []; spY = []; spAlloc = [];
for i = 1:nc-1
    % extract flow field
    singleField = squeeze(v(:,:,i,:));
    
    % create warping operator forward
    warp = warpingOperator([Ny,Nx],singleField);
    idOp = speye(Nx*Ny);
    
    % find out of range warps in each of the operators and set the corresponding line in the other operator also to zero
    marker = sum(abs(warp),2) == 0;
    
    % find occlusion markers from optional occlusion map:
    if ~isempty(occ)
        marker_occ = sparse(vec(occ(:,:,i)));
        marker = logical(marker+marker_occ); 
    end
    
    warp(marker > 0,:) = 0;
    idOp(marker > 0,:) = 0; %#ok<SPRIX>
    
    if strcmp(type,'I-F')
        spX = [spX;((Nx*Ny*(i-1)+1):(Nx*Ny*i))']; %#ok<*AGROW>
        spY = [spY;((Nx*Ny*(i-1)+1):(Nx*Ny*i))'];
        spAlloc = [spAlloc;-full(diag(idOp))];
        
        [sp_x,sp_y,walloc] = find(warp);
        spX = [spX;(Nx*Ny*(i-1)+1)+sp_x-1];
        spY = [spY;(Nx*Ny*i+1)+sp_y-1];
        spAlloc = [spAlloc;walloc];
    elseif strcmp(type,'F-I')
        [sp_x,sp_y,walloc] = find(warp);
        spX = [spX;(Nx*Ny*(i-1)+1)+sp_x-1];
        spY = [spY;(Nx*Ny*(i-1)+1)+sp_y-1];
        spAlloc = [spAlloc;walloc];
        
        spX = [spX;((Nx*Ny*(i-1)+1):(Nx*Ny*i))'];
        spY = [spY;((Nx*Ny*i+1):(Nx*Ny*(i+1)))'];
        spAlloc = [spAlloc;-full(diag(idOp))];
    end
    
end
warpingOp = sparse(spX,spY,spAlloc,Nx*Ny*nc,Nx*Ny*nc);
end


function x = vec(x)
    x = x(:);
end
