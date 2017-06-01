%% MULTIFRAME MOTION COUPLING FOR VIDEO SUPER RESOLUTION
% 
% 
% Be sure to initialize the submodules and to compile them
% following their instructions


clearvars;
%addpath(genpath(cd)); 



%% Data properties

startFrame = 1;
numFrames = 18; 
cslice = ceil(numFrames/2);

factor  = 4;             % Magnification factor


%% Load VSR data
% city:
datasetName = '/windows/DataJonas/ScieboLocalFolder/Data/VSR_Sun_dropbox/vsr_result/city_org_I420/input'; % 33 frames
gtName = '/windows/DataJonas/ScieboLocalFolder/Data/VSR_Sun_dropbox/vsr_result/city_org_I420/original';

% walk:
%datasetName = '/windows/DataJonas/ScieboLocalFolder/Data/VSR_Sun_dropbox/vsr_result/people/input'; % 47 frames
%gtName = '/windows/DataJonas/ScieboLocalFolder/Data/VSR_Sun_dropbox/vsr_result/people/original';

% calendar:
%datasetName = '/windows/DataJonas/ScieboLocalFolder/Data/VSR_Sun_dropbox/vsr_result/mobile/input';
%gtName = '/windows/DataJonas/ScieboLocalFolder/Data/VSR_Sun_dropbox/vsr_result/mobile/original';


% Load
[imageSequenceSmall] = LoadImSequenceReal(datasetName,startFrame,numFrames);   
[~,imageSequenceLarge] = LoadImSequence(gtName,startFrame,numFrames,factor,'bicubic');   




%% Construct algorithm object
% Input: RGB-Time matlab array 
mainSuper = MultiframeMotionCoupling(imageSequenceSmall);

%% Set variables (these are the standard parameters)

% Procedure
mainSuper.factor        = factor;              % magnification factor
mainSuper.verbose       = 1;                   % enable intermediate output, 1 is text, 2 is image
mainSuper.framework     = 'prost';           % Choose framework for super resolution problem
                                               % Either 'flexBox' or 'flexBox_vector'
                                               % or 'prost', if installed

% Problem parameters
mainSuper.alpha         = 0.001;                % regularizer weight
mainSuper.beta          = 0.2;                 % flow field complexity
mainSuper.kappa         = 0.1;                % regularization pendulum
mainSuper.flowDirection = 'forward';           % flow field direction

% Operator details
mainSuper.interpMethod = 'average';            % Downsampling operator D'average';%
mainSuper.k = fspecial('gaussian',7,sqrt(0.6));% Blur operator B  sqrt(0.6));%
mainSuper.VDSRFlag     = false;                 % Enable VDSR initial guess


%% Init flow field and solvers

if mainSuper.VDSRFlag
    vl_setupnn();
end
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
prost.release

%% Check against LiuSun psnr:
%
% check LiuSun PSNR
% 26.03 for frame 28


%% Load VSR data
% city:
outputName = '/windows/DataJonas/ScieboLocalFolder/Data/VSR_Sun_dropbox/vsr_result/city_org_I420/output'; % 33 frames

% Load
[~,imageSequenceVSR] = LoadImSequence(outputName,startFrame,numFrames,factor,'bicubic');
disp_('VSR PSNR :')
for ii = 1:numFrames
    outImage = imageSequenceVSR(20:end-20,20:end-20,:,ii);
    psnrErr(ii) = round(psnr(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ii)),2);
    ssimErr(ii) = round(ssim(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ii)),3);
end
figure(1),plot(startFrame:startFrame+numFrames-1,psnrErr), title(['PSNR (mean): ',num2str(mean(psnrErr))]);
figure(2),plot(startFrame:startFrame+numFrames-1,ssimErr),title(['SSIM (mean):',num2str(mean(ssimErr))]);