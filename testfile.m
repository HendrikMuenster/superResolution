%% MULTIFRAME MOTION COUPLING FOR VIDEO SUPER RESOLUTION
% 
% 
% Be sure to initialize the submodules and to compile them
% following their instructions


clearvars;
%addpath(genpath(cd)); 



%% Data properties
datasetName = 'city';

startFrame = 1;
numFrames = 13;
cslice = ceil(numFrames/2);


factor = 4;       % Magnification factor

%% ICCV paper data generation process 
dataFolder = '/windows/DataJonas/ScieboLocalFolder/Data/videos_scenes/';
[imageSequenceSmall,imageSequenceLarge] = LoadImSequence([dataFolder,filesep,datasetName],startFrame,numFrames,factor,'bicubic');   



%% Construct algorithm object
% Input: RGB-Time matlab array 
mainSuper = MultiframeMotionCoupling(imageSequenceSmall);

%% Set variables (these are the standard parameters)

% Procedure
mainSuper.factor        = factor;              % magnification factor
mainSuper.verbose       = 1;                   % enable intermediate output, 1 is text, 2 is image

% Problem parameters
mainSuper.alpha         = 0.01;                % regularizer weight
mainSuper.beta          = 0.2;                 % flow field complexity
mainSuper.kappa         = 0.25;                % regularization pendulum
mainSuper.flowDirection = 'forward';           % flow field direction

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
     vid = implay(mainSuper.result1,2);  
     set(vid.Parent, 'Position',get(0, 'Screensize'));
else
    figure(1), imshow(outImage); title(['PSNR: ', num2str(psnrErr)]); axis image
end
%% 
disp('---------------------------------------------------------------------')
