function [imgSR,timings] = singleframeMotionSR_unger_prost(imageSequenceSmall,factor,alpha,beta,huber_eps,comp_mode)
% Unger - Werlberger algorithm
% as described in
% Unger, Markus, et al. "A Convex Approach for Variational Super-Resolution." DAGM-Symposium. 2010.
%



%% 1) Init
timings = zeros(2,1);
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
    uTmp = cat(3,J1(:,:,i),J1(:,:,central_slice));
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
    v(:,:,i,2) = vTmp(:,:,1,1);
    v(:,:,i,1) = vTmp(:,:,1,2);
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
backend = prost.backend.pdhg('stepsize', 'boyd','tau0', 1, ...
    'sigma0', 1);                       % prost backend options
%obj.opts.backend = prost.backend.pdhg('stepsize', 'goldstein', ...
%                 'residual_iter', 100);
if strcmp(comp_mode,'accurate')
    iters = 7500;
else
    iters = 2500;
end
opts = prost.options('max_iters', iters, 'num_cback_calls', 5,...
    'verbose', true);

% Add primal variable:
N = yLarge*xLarge;
n = ySmall*xSmall;
u = prost.variable(N);
% Add dual variables
p = prost.variable(n*numFrames);

if numFrames == 5 % whyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
    p1 = prost.sub_variable(p,n);
    p2 = prost.sub_variable(p,n);
    p3 = prost.sub_variable(p,n);
    p4 = prost.sub_variable(p,n);
    p5 = prost.sub_variable(p,n);
elseif numFrames == 13
    p1 = prost.sub_variable(p,n);
    p2 = prost.sub_variable(p,n);
    p3 = prost.sub_variable(p,n);
    p4 = prost.sub_variable(p,n);
    p5 = prost.sub_variable(p,n);
    p6 = prost.sub_variable(p,n);
    p7 = prost.sub_variable(p,n);
    p8 = prost.sub_variable(p,n);
    p9 = prost.sub_variable(p,n);
    p10= prost.sub_variable(p,n);
    p11= prost.sub_variable(p,n);
    p12= prost.sub_variable(p,n);
    p13= prost.sub_variable(p,n);
else
    error('todo');
end
q = prost.variable(2*N); % TV!
SMC = prost.min_max_problem({u},{p,q});

% Add data term pairings
SMC.add_dual_pair(u,p1,prost.block.sparse(D*B*warp{1}));
SMC.add_dual_pair(u,p2,prost.block.sparse(D*B*warp{2}));
SMC.add_dual_pair(u,p3,prost.block.sparse(D*B*warp{3}));
SMC.add_dual_pair(u,p4,prost.block.sparse(D*B*warp{4}));
SMC.add_dual_pair(u,p5,prost.block.sparse(D*B*warp{5}));
if numFrames == 13 % this code feels violating ;< 
    SMC.add_dual_pair(u,p6,prost.block.sparse(D*B*warp{6}));
    SMC.add_dual_pair(u,p7,prost.block.sparse(D*B*warp{7}));
    SMC.add_dual_pair(u,p8,prost.block.sparse(D*B*warp{8}));
    SMC.add_dual_pair(u,p9,prost.block.sparse(D*B*warp{9}));
    SMC.add_dual_pair(u,p10,prost.block.sparse(D*B*warp{10}));
    SMC.add_dual_pair(u,p11,prost.block.sparse(D*B*warp{11}));
    SMC.add_dual_pair(u,p12,prost.block.sparse(D*B*warp{12}));
    SMC.add_dual_pair(u,p13,prost.block.sparse(D*B*warp{13}));
end

% Add regularizer
SMC.add_dual_pair(u,q,prost.block.gradient2d(xLarge,yLarge,1));

% Add functions
SMC.add_function(p1, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, vec(J0(:,:,1)), huber_eps));
SMC.add_function(p2, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, vec(J0(:,:,2)), huber_eps));
SMC.add_function(p3, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, vec(J0(:,:,3)), huber_eps));
SMC.add_function(p4, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, vec(J0(:,:,4)), huber_eps));
SMC.add_function(p5, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, vec(J0(:,:,5)), huber_eps));
if numFrames == 13 %uuugh
    SMC.add_function(p6, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, vec(J0(:,:,6)), huber_eps));
    SMC.add_function(p7, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, vec(J0(:,:,7)), huber_eps));
    SMC.add_function(p8, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, vec(J0(:,:,8)), huber_eps));
    SMC.add_function(p9, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, vec(J0(:,:,9)), huber_eps));
    SMC.add_function(p10, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, vec(J0(:,:,10)), huber_eps));
    SMC.add_function(p11, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, vec(J0(:,:,11)), huber_eps));
    SMC.add_function(p12, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, vec(J0(:,:,12)), huber_eps));
    SMC.add_function(p13, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, vec(J0(:,:,13)), huber_eps));
end
SMC.add_function(q, prost.function.sum_norm2(2, false, 'ind_leq0', 1/alpha, 1, 1, 0, huber_eps/alpha)); %l^{2,1}

disp('Single frame motion coupling initialized')
%tic;
tic
prost.solve(SMC,backend,opts);
timings(2) = toc;
%toc;

I = reshape(u.val,yLarge,xLarge);

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