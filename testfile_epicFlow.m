%% MULTIFRAME MOTION COUPLING FOR VIDEO SUPER RESOLUTION
% 
% 
% Be sure to initialize the submodules and to compile them
% following their instructions
%
% 
% This is an example file for the incorporation of alternative
% optical flow algorithms, here a MATLAB wrapper for epicFlow
% [ https://github.com/suhangpro/epicflow ]




clearvars;
addpath(genpath(cd)); 



%% Data properties
datasetName = 'surfer';

startFrame = 5;
numFrames = 5;
cslice = ceil(numFrames/2);


factor = 4;       % Magnification factor

%% Data generation process 
dataFolder = '/windows/DataJonas/ScieboLocalFolder/Data/videos_scenes/';
[imageSequenceSmall,imageSequenceLarge] = LoadImSequence([dataFolder,filesep,datasetName],startFrame,numFrames,factor,'bicubic');   


%% Compute epicFlow
% epicFlow parameters:
%   interpolation parameters
%     -nw                                                      use Nadaraya-Watson instead of LA interpolator in the interpolation
%     -p, -prefnn             <int>(25)                        number of neighbors for consisteny checking in the interpolation
%     -n, -nn                 <int>(100)                       number of neighnors for the interpolation
%     -k                      <float>(0.8)                     coefficient of the sigmoid of the Gaussian kernel used in the interpolation
%   energy minimization parameters
%     -i, -iter               <int>(5)                         number of iterations for the energy minimization
%     -a, -alpha              <float>(1.0)                     weight of smoothness term
%     -g, -gamma              <float>(3.0)                     weight of gradient constancy assumption
%     -d, -delta              <float>(2.0)                     weight of color constancy assumption
%     -s, -sigma              <float>(0.8)                     standard deviation of Gaussian presmoothing kernel
%   predefined parameters
%     -sintel                                                  set the parameters to the one optimized on (a subset of) the MPI-Sintel dataset
%     -middlebury                                              set the parameters to the one optimized on the Middlebury dataset
%     -kitti                                                   set the parameters to the one optimized on the KITTI dataset
params = '-prefnn 250 -k 0.6 -i 25 -a 2 -s 0.6';
    
v = zeros(ny*factor,nx*factor,nf-1,2);
for i = 1:numFrames-1
    if mod(i,2) == 1
        get_epicflow(imageSequenceSmall(:,:,:,i),imageSequenceSmall(:,:,:,i+1),'temp_flow.flo',params);
        vTmp = readFlowFile('temp_flow.flo');
        v(:,:,i,1) = imresize(vTmp(:,:,1),factor) * factor;
        v(:,:,i,2) = imresize(vTmp(:,:,2),factor) * factor;
    else
        get_epicflow(imageSequenceSmall(:,:,:,i+1),imageSequenceSmall(:,:,:,i),'temp_flow.flo',params);
        vTmp = readFlowFile('temp_flow.flo');
        v(:,:,i,1) = imresize(vTmp(:,:,1),factor) * factor;
        v(:,:,i,2) = imresize(vTmp(:,:,2),factor) * factor; 
    end
    figure(i),imagesc(flowToColorV2(vTmp)), drawnow
end


%% Construct algorithm object
% Input: RGB-Time matlab array 
mainSuper = MultiframeMotionCoupling(imageSequenceSmall,'flowField',v);

%% Set variables

% Procedure
mainSuper.factor        = factor;              % magnification factor
mainSuper.verbose       = 1;                   % enable intermediate output, 1 is text, 2 is image

% Problem parameters
mainSuper.alpha         = 0.01;                % regularizer weight
mainSuper.beta          = 0.2;                 % flow field complexity
mainSuper.kappa         = 0.25;                % regularization pendulum
mainSuper.flowDirection = 'backward';           % flow field direction

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
