%% MULTIFRAME MOTION COUPLING FOR VIDEO SUPER RESOLUTION
% 
% 
% Be sure to initialize the submodules and to compile them
% following their instructions
% Total equivalence of forward and backward model, with an easy downsampling procedure.
% 
%% Data properties

startFrame = 23;
numFrames = 13; 
cslice = ceil(numFrames/2);

factor  = 4;             % Magnification factor


%% Load VSR data
dataLoc = '../data/video_scenes/city';

  

% Generate data with 
interpMethod = 'stride';            % Downsampling operator D
k = fspecial('gaussian',7,1.6);       % Blur operator B  
Quantize = false;
[imageSequenceSmall,imageSequenceLarge] = LoadImSequenceForward(dataLoc,startFrame,numFrames,factor,interpMethod,k,Quantize);  


%% Construct algorithm object
% Input: RGB-Time matlab array 
mainSuper = MultiframeMotionCoupling(imageSequenceSmall);

%% Set variables (these are the standard parameters)

% Procedure
mainSuper.factor        = factor;              % magnification factor
mainSuper.verbose       = 1;                   % enable intermediate output, 1 is text, 2 is image
mainSuper.framework     = 'prost';             % Choose framework for super resolution problem
                                               % Either 'flexBox' or 'flexBox_vector'
                                               % or 'prost', if installed

% Problem parameters
mainSuper.alpha         = 0.05;                % regularizer weight
mainSuper.beta          = 1;                 % flow field complexity
mainSuper.kappa         = NaN;                 % regularization balance
mainSuper.flowDirection = 'backward';          % flow field direction

% Operator details
mainSuper.interpMethod = interpMethod;         % Downsampling operator D'average';%
mainSuper.k = k;                               % Blur operator B 



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
     vid = implay(mainSuper.result1);  
     set(vid.Parent, 'Position',get(0, 'Screensize'));
else
    figure(1), imshow(outImage); title(['PSNR: ', num2str(psnrErr)]); axis image
end
%% 
disp('---------------------------------------------------------------------')

