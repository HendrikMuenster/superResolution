%
% Video Super Resolution

clearvars;


%% Data properties
datasetName = 'city';

startFrame = 1;
numFrames = 13;
cslice = ceil(numFrames/2);


factor = 4;                         % magnification factor

%% Load Video and code
dataFolder = '/windows/DataJonas/ScieboLocalFolder/Data/videos_scenes/';
[imageSequenceSmall,imageSequenceLarge] = LoadImSequence([dataFolder,filesep,datasetName],startFrame,numFrames,factor,'bicubic',0);   



%% Run Single frame motion coupling
alpha = 0.05;
beta = 0.1;

imgSR = singleframeMotionSR(imageSequenceSmall,factor,alpha,beta);

%% Show error margin

outImage = imgSR(20:end-20,20:end-20,:);
psnrErr = psnr(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2)));
ssimErr = ssim(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2)));
disp(['PSNR : ',num2str(psnrErr),' dB']);
disp(['SSIM : ',num2str(ssimErr),' ']);

figure(1), imshow(outImage); title(['PSNR: ', num2str(psnrErr)]); axis image
disp('---------------------------------------------------------------------')


%% Compare total psnr/ssim

