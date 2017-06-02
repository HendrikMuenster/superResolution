%% MULTIFRAME MOTION COUPLING FOR VIDEO SUPER RESOLUTION
% 
% 
% Be sure to initialize the submodules and to compile them
% following their instructions


clearvars;
%addpath(genpath(cd)); 
cd ~/MATLAB/superResolution


%% Data properties

startFrame = 15;
numFrames = 18; 
cslice = ceil(numFrames/2);

factor  = 4;             % Magnification factor


%% Load VSR data
gtName = '/windows/DataJonas/ScieboLocalFolder/Data/VSR_Sun_dropbox/vsr_result/city_org_I420/original';
%gtName = '/windows/DataJonas/ScieboLocalFolder/Data/VSR_Sun_dropbox/vsr_result/people/original';
%gtName = '/windows/DataJonas/ScieboLocalFolder/Data/VSR_Sun_dropbox/vsr_result/mobile/original';

datasetName = '/windows/DataJonas/ScieboLocalFolder/Data/VSR_Sun_dropbox/vsr_result/city_org_I420/input'; % 33 frames
[imageSequenceSmall] = LoadImSequenceReal(datasetName,startFrame,numFrames);   
% Load forward results
interpMethod = 'stride';            % Downsampling operator D
k = fspecial('gaussian',7,1.6);       % Blur operator B  
Quantize = true;
%[imageSequenceSmall,imageSequenceLarge] = LoadImSequenceForward(gtName,startFrame,numFrames,factor,interpMethod,k,Quantize);  


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
mainSuper.alpha         = 0.01;                % regularizer weight
mainSuper.beta          = 0.2;                 % flow field complexity
mainSuper.kappa         = 0.25;                   % regularization pendulum
mainSuper.flowDirection = 'backward';          % flow field direction

% Operator details
mainSuper.interpMethod = interpMethod;         % Downsampling operator D'average';%
mainSuper.k = k;                               % Blur operator B 
mainSuper.VDSRFlag     = false;                % Enable VDSR initial guess


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
     vid = implay(mainSuper.result1);  
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
% 26.03 for frame 23


%% Compare total psnr/ssim
for ii = 1:numFrames
    outImage = mainSuper.result1(20:end-20,20:end-20,:,ii);
    psnrErr(ii) = round(psnr(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ii)),2);
    %ssimErr(ii) = round(ssim(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ii)),3);
end
uiopen('/home/vsa_jonas/MATLAB/superResolution/liuSun_psnr_city.fig',1), hold on
plot(startFrame:startFrame+numFrames-1,psnrErr), title(['PSNR (mean): ',num2str(mean(psnrErr))]);
%figure(2),plot(startFrame:startFrame+numFrames-1,ssimErr),title(['SSIM (mean):',num2str(mean(ssimErr))]);
legend('LiuSun','Ours','Location','northwest')
hold off
