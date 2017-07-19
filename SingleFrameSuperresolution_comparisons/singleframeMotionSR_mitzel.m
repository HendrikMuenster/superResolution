function imgSR = singleframeMotionSR_mitzel(imageSequenceSmall,factor,alpha,beta)
% Unger - Werlberger algorithm
%



%% 1) Init
%tic
[ySmall,xSmall,nColors,numFrames] = size(imageSequenceSmall);
yLarge = ySmall*factor; xLarge = xSmall*factor;

% Convert to YCbCr
if nColors > 1
    J0 = zeros(ySmall,xSmall,numFrames);
    for i = 1:numFrames
        YcbCr = rgb2ycbcr(squeeze(imageSequenceSmall(:,:,:,i)));
        J0(:,:,i) = YcbCr(:,:,1);
    end
end
central_slice = ceil(numFrames/2);

% Get bicubic upsamplings
J1 = zeros(yLarge,xLarge,numFrames);
for i = 1:numFrames
    J1(:,:,i) = imresize(J0(:,:,i),factor,'bicubic');
end

% Get initial flow fields from motionEstimator
v = zeros(yLarge,xLarge,numFrames,2);
disp('Flow precomputation initialized...');

for i = 1:numFrames
    if i == central_slice
        continue
    end
    uTmp = cat(3,J0(:,:,i),J0(:,:,central_slice));
    motionEstimator = motionEstimatorClass(uTmp,1e-6,beta);
    % optical flow parameters
    motionEstimator.verbose = 0;
    motionEstimator.dataTerm = 'L1';
    motionEstimator.regularizerTerm = 'Huber';
    motionEstimator.doGradientConstancy = 1;
    motionEstimator.doGaussianSmoothing = 1;
    motionEstimator.medianFiltering = 1;
    motionEstimator.steplength = 0.8;
    motionEstimator.numberOfWarps = 3;
    % run OF
    motionEstimator.init;
    motionEstimator.runPyramid;
    
    vTmp = motionEstimator.getResult;
    v(:,:,i,2) = factor*imresize(vTmp(:,:,1,1),factor);
    v(:,:,i,1) = factor*imresize(vTmp(:,:,1,2),factor);
    disp(['Optical Flow to frame ',num2str(i),' computed']);
end

% Compute warp operators

for i = 1:numFrames
    if i == central_slice
        warp{i} = speye(yLarge*xLarge); %#ok<*AGROW>
    end
    warp{i} = warpingOperator([yLarge,xLarge],squeeze(v(:,:,i,:))); 
end



% Initialize sampling and blur matrices
width = 7;
kernel = @(x) double(abs(x)<=2)/4;
D = samplingOperator([yLarge,xLarge],[ySmall,xSmall],'custom',{kernel,width},false);
sigm_Unger = 0.25*sqrt(factor^2-1);
B  = RepConvMtx(fspecial('gaussian',7,sigm_Unger),[yLarge,xLarge]);
%toc 

%% 2) Build flexBox object and solve
%tic
SMC = flexBox;
SMC.params.tol = 1e-4;
SMC.params.tryCPP = 1;
SMC.params.verbose = 1;
SMC.params.maxIt = 10000;
h_eps = 0.01;

% Add primal variable:
u_id = SMC.addPrimalVar([yLarge,xLarge]);

% Add data terms
for i = 1:numFrames
    SMC.addTerm(L1dataTermOperator(1,D*B*warp{i},J0(:,:,i)),u_id);
end

% Add regularizer
SMC.addTerm(L1gradientIso(alpha,[yLarge,xLarge]),u_id);%

%Grad = spmat_gradient2d(yLarge,xLarge,1);
%SMC.addTerm(huberDataTermOperator(alpha,Grad,zeros(yLarge,xLarge,2),h_eps),u_id);

disp('Single frame motion coupling initialized')
%tic;
SMC.runAlgorithm;
%toc;

I = SMC.getPrimal(u_id);

%toc
%% 3) Finalize

% Reconvert to YCbCr
I_ycbcr = imresize(rgb2ycbcr(imageSequenceSmall(:,:,:,central_slice)),factor,'bicubic');
I_ycbcr(:,:,1) = I;



% Return output I_sr

imgSR = ycbcr2rgb(I_ycbcr);
end




%%%%% subfunctions

function u_vec  = vec(u)
% vectorize multidimensional u

u_vec = u(:);
end