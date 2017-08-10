%
% Video Super Resolution

clearvars;


%% Data properties
datasetName = 'city';

startFrame = 1;
numFrames = 3;
cslice = ceil(numFrames/2);


factor = 4;                         % magnification factor

%% Load Video and code
dataFolder = 'data/video_scenes/';
[imageSequenceSmall,imageSequenceLarge] = LoadImSequence([dataFolder,filesep,datasetName],startFrame,numFrames,factor,'bicubic',0);   



%% Run Single frame motion coupling
alpha = 0.1;
beta = 0.2;
tic
imgSR = singleframeMotionSR_unger_prost(imageSequenceSmall,factor,alpha,beta,'accurate');
toc
%% Show error margin

outImage = imgSR(20:end-20,20:end-20,:);
psnrErr = psnr(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2)));
ssimErr = ssim(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2)));
disp(['PSNR : ',num2str(psnrErr),' dB']);
disp(['SSIM : ',num2str(ssimErr),' ']);

figure(1), imshow(outImage); title(['PSNR: ', num2str(psnrErr)]); axis image
disp('---------------------------------------------------------------------')
