function [S] = samplingOperator(inputDims,targetDims,method,kernel,antialiasing)
%
% Compute a sampling operator that samples a 2D signal with given inputDims
% to targetDims.
%
% The type of interpolation is given by 'method'.
% Valid standard methods are 'nearest', 'bilinear','bicubic',lanczos3','pchip'.
%
% If the method is 'custom', then a function handle to a kernel can be
% given in kernel, written as {function,width} in analogy to the imresize function
% by for example {@lanczos2,4.0} with
% % function f = lanczos2(x)
% % f = (sin(pi*x) .* sin(pi*x/2) + eps) ./ ((pi^2 * x.^2 / 2) + eps);
% % f = f .* (abs(x) < 2);
% % end
%
%
% 'nearest'     -> nearest neighbors
% 'bilinear'    -> linear interpolation
% 'bicubic'     -> cubic spline interpolation, based on imresize
% 'bicubic-0'   -> cubic spline interpolation, based on interp1 (faster construction but zero-padding)
% 'lanczos2'    -> lanczos kernel based interpolation a = 2
% 'lanczos3'    -> lanczos kernel based interpolation a = 3
% 'pchip'       -> shape preserving piecewise cubic Hermite interpolating polynomial ('local cubic')
% 'average'     -> 'averaging' block filter, taking the average over all pixels in the superpixel
% 'stride'      -> the 'machine learning thing'. Start in the upper left corner and take each factor-th pixel
%
%
% The boolean antialiasing calls the MATLAB antialiasing routine which widens the interpolation kernel
% custom anti-aliasing kernels can be included in a custom kernel.
%
% Implementation roughly follows Michael's lecture note after the
% discovery that the MATLAB GIT compiler makes this approach not as
% terribly inefficient one would expected for imresize and can be
% parallelized for interp1 calls
% JOG


% Validate interpolation method
validatestring(method,{'nearest','bilinear','bicubic','bicubic-0','lanczos2','lanczos3','pchip','custom','average','stride'});


% Validate input arguments
if nargin < 5
    antialiasing = false;
end

% Bilinear and fast bicubic will be called differently in interp1 code
if strcmp(method,'bilinear') && ~antialiasing
    method = 'linear';
elseif strcmp(method,'bicubic-0') && ~antialiasing
    method = 'v5cubic';
end


if nargin < 4
    if strcmp(method,'custom')
        error('No custom kernel given');
    else
        kernel = [];
    end
end

% pchip and fast bicubic are not possible with antialiasing
if antialiasing && ismember(method,{'v5cubic','pchip'})
    error(' The chosen method is not possible together with antialiasing');
end

% average is shorthand for this custom filter:
if strcmp(method,'average')
    method = 'custom';
    kernel = {@(x) double(abs(x)<=2)/4,7}; % this is the average kernel with 7 pixels
end



% Compute sampling matrices Sx and Sy separately:
Sx = Dim1Matrix(inputDims(2),targetDims(2),method,kernel,antialiasing) ;
Sy = Dim1Matrix(inputDims(1),targetDims(1),method,kernel,antialiasing)';


% Call kronecker product to compute output matrix
S = kron(Sx',Sy);

end

%=====================================================================

function S = Dim1Matrix(inputDim,targetDim,method,kernel,antialiasing)
% Compute actual 1D sampling matrix
% Input as above
%

% Compute sampling factor
factor = inputDim/targetDim;


% Call either interp1 (fast) or imresize(slow) to build kernel
if ismember(method,{'nearest','linear','pchip','v5cubic'}) && ~antialiasing
    
    % Construct input and output sampling positions
    inputPos = (1:inputDim);
    targetPos = (1+factor/2:factor:inputDim+1) - 0.5;
    
    % Call interp1 with the identity matrix
    S = sparse(interp1(inputPos,eye(inputDim),targetPos,method))';
    
elseif ismember(method,{'nearest','bilinear','bicubic','lanczos2','lanczos3'}) || antialiasing && ~strcmp(method,'custom')
    
    % Allocate matrix with enough memory
    if antialiasing
        S = spalloc(inputDim,targetDim,4*factor*targetDim);
    else
        S = spalloc(inputDim,targetDim,2*factor*targetDim);
    end
    
    % Iterate over all dimensions and call imresize separately
    for ii = 1:inputDim % this loop is bearable thanks to the GIT compiler
        eyeVector = zeros(1,inputDim); eyeVector(ii) = 1; % ii-th row of 2D identity
        S(ii,:) = imresize(eyeVector,1/factor,method,'Antialiasing',antialiasing); %#ok<SPRIX>
    end
    
elseif strcmp(method,'custom')
    
    % Allocate matrix with enough memory
    if antialiasing
        S = spalloc(inputDim,targetDim,4*factor*targetDim);
    else
        S = spalloc(inputDim,targetDim,2*factor*targetDim);
    end
    
    % Iterate over all dimensions and call imresize separately
    for ii = 1:inputDim % this loop is bearable thanks to the GIT compiler
        eyeVector = zeros(1,inputDim); eyeVector(ii) = 1; % ii-th row of 2D identity
        S(ii,:) = imresize(eyeVector,1/factor,kernel,'Antialiasing',antialiasing); %#ok<SPRIX>
    end
elseif strcmp(method,'stride')
    
    mcords = 1:factor:inputDim;
    ncords = 1:targetDim;
    S      = sparse(mcords,ncords,1,inputDim,targetDim); % is actually S'
    
    
end
%disp(S(:,87)); % for debugging
end

