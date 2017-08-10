%% MULTIFRAME MOTION COUPLING FOR VIDEO SUPER RESOLUTION
% 
% 
% Be sure to initialize the submodules and to compile them
% following their instructions


clearvars;



%% Data properties
dataFolder = 'data/video_scenes/';
flowFolder = 'data/flow_scenes/';
occFolder = 'data/occlusion_scenes/';
datasetName = 'bandage_1_final';
datasetName_flow = 'bandage_1';

startFrame = 1;
numFrames = 13;

factor  = 4;             % Magnification factor



%% Data generation process 

[imageSequenceSmall,imageSequenceLarge] = LoadImSequence([dataFolder,filesep,datasetName],startFrame,numFrames,factor,'bicubic');   

%% load gt flow


[ny,nx,nc,nf] = size(imageSequenceSmall);
v = zeros(ny*factor,nx*factor,nf-1,2);
for i = 1:nf-1
    vTmp = readFlowFile([flowFolder,filesep,datasetName_flow,filesep,'frame_',num2str(i,'%04d'),'.flo']);
    v(:,:,i,1) = vTmp(:,:,1);
    v(:,:,i,2) = vTmp(:,:,2);
end

%% load occclusions
occ = zeros(ny*factor,nx*factor,nf-1);
for i = 1:nf-1
    occTmp = im2double(imread([occFolder,filesep,datasetName_flow,filesep,'frame_',num2str(i,'%04d'),'.png']));
    occ(:,:,i) = logical(occTmp);
end

%% construct occluded warp operator
flowDir = 'backward';

% Construct warp in specified direction with occlusions
if strcmp(flowDir,'forward')
    warpingOp = constructWarpFF(v,'F-I',occ);
elseif strcmp(flowDir,'backward')
    warpingOp = constructWarpFF(v,'I-F',occ);
elseif strcmp(flowDir,'forward-backward')
    warpingOp = constructWarpFB(v,occ);
end


%% Construct algorithm object
% Input: RGB-Time matlab array 

%mainSuper = MultiframeMotionCoupling(imageSequenceSmall);               % standard
mainSuper = MultiframeMotionCoupling(imageSequenceSmall,'warpOp',warpingOp); % GT flow with occlusions known
%mainSuper = MultiframeMotionCoupling(imageSequenceSmall,'flowField',v);% GT flow without occlusions
%% Set variables (these are the standard parameters)

% Procedure
mainSuper.factor        = factor;              % magnification factor
mainSuper.verbose       = 1;                   % enable intermediate output, 1 is text, 2 is image
mainSuper.framework     = 'prost';             % Choose framework for super resolution problem
                                               % Either 'flexBox' or 'flexBox_vector'
                                               % or 'prost', if installed

% Problem parameters
mainSuper.alpha         = 0.01;                % regularizer weight
mainSuper.beta          = 0.2;                 % flow field complexity
mainSuper.kappa         = 0.5;                 % regularization pendulum
mainSuper.flowDirection = flowDir;             % flow field direction 

% Operator details
mainSuper.interpMethod = 'average';            % Downsampling operator D
mainSuper.k = fspecial('gaussian',7,sqrt(0.6));% Blur operator B  


%% Init flow field and solvers
tic
mainSuper.init;
toc

%% Solve super resolution problem
tic
mainSuper.run;
toc

%% Show error margin

outImage = mainSuper.result1(20:end-20,20:end-20,:,ceil(numFrames/2));
psnrErr = round(psnr(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2))),2);
ssimErr = round(ssim(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2))),3);
disp(['PSNR (central patch, central slice): ',num2str(psnrErr),' dB']);
disp(['SSIM (central patch, central slice): ',num2str(ssimErr),' ']);


%% Visualize either central image or full video
if mainSuper.verbose > 0
     %vid = implay(mainSuper.result1,2);  
     %set(vid.Parent, 'Position',get(0, 'Screensize'));
else
    figure(1), imshow(outImage); title(['PSNR: ', num2str(psnrErr)]); axis image
end