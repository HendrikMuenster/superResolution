%
% Video Super Resolution

clearvars;

%% Data properties
datasetName = 'penguins';

startFrame = 1;
numFrames = 5;
cslice = ceil(numFrames/2);


factor = 4;                         % magnification factor

%% Load Video and code
dataFolder = '/windows/DataJonas/ScieboLocalFolder/Data/videos_scenes/';
[imageSequenceSmall,imageSequenceLarge] = LoadImSequence([dataFolder,filesep,datasetName],startFrame,numFrames,factor);    
addpath(genpath(cd)); % only load paths if data location was error free
[ny,nx,nc,nf] = size(imageSequenceSmall);

%% Compute epicFlow and load into MMC
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

%% Init algorithm class thing
mainSuper = MultiframeMotionCoupling(imageSequenceSmall,'flowField',v);

%% Set variables

% Procedure
mainSuper.factor        = factor;              % magnification factor
mainSuper.numMainIt     = 1;                   % number of total outer iterations
mainSuper.verbose       = 1;                   % enable intermediate output, 1 is text, 2 is image

% Problem parameters
mainSuper.alpha1        = 0.01;                % temporal weight
mainSuper.alpha2        = 0.01;                % spatial weight
mainSuper.beta          = 0.1;                 % flow field complexity
mainSuper.kappa         = 0.2;                 % regularization pendulum

  
% Downsampling details
mainSuper.interpMethod = 'custom';
mainSuper.interpAA     = false;                % Increases the scaling of the kernel to prevent aliasing - MATLAB inverse, handle with care(!)

if strcmp(mainSuper.interpMethod,'custom')
    width = 7;
    kernel = @(x) double(abs(x)<=2)/4; % this is the old average kernel
    mainSuper.interpKernel = {kernel,width};
end

% "Motion" blur:
mainSuper.k = fspecial('gaussian',7,sqrt(0.6));

% %% Plot kernel comparisons
% h = 0.0001;
% intervalThing = -width:h:width;
% sigSTD = sqrt(1.2/2);
% gaussThing = @(x) 1./sqrt(2*pi*sigSTD^2).*exp(-x.^2 / 2/sigSTD^2);
% figure(1), plot(intervalThing,conv(gaussThing(intervalThing),kernel(intervalThing),'same')*h);
% hold on, plot(intervalThing,gaussThing(intervalThing))
% plot(intervalThing,kernel(intervalThing)),
% 
% sigSTD = 1.6;
% gaussThing = @(x) 1./sqrt(2*pi*sigSTD^2).*exp(-x.^2 / 2/sigSTD^2);
% plot(intervalThing,gaussThing(intervalThing));
% legend('Our kernel','Gaussian sigma=0.776 ','Block sampling','Gaussian sigma=1.6 ')


%% INIT flow field and solvers
mainSuper.init;


%% Solve super resolution problem

mainSuper.run;


%% Show error margin

outImage = mainSuper.result1(20:end-20,20:end-20,:,ceil(numFrames/2));
psnrErr = psnr(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2)));
ssimErr = ssim(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2)));
disp(['PSNR (central patch, central slice): ',num2str(psnrErr),' dB']);
disp(['SSIM (central patch, central slice): ',num2str(ssimErr),' ']);

%% visualize video
if mainSuper.verbose > 0
    vid = implay(mainSuper.result1,2);  % persists through close all ( \_('')_/� ) - use this to compare visual quality to previous iterations
    set(vid.Parent, 'Name', ['u, computed with kappa = ',num2str(mainSuper.kappa),...
                            ' and alphas ',num2str(mainSuper.alpha1),', ',num2str(mainSuper.alpha2),...
                            ', OF is Huber-L1 with beta = ',num2str(mainSuper.beta)]);
    set(vid.Parent, 'Position',get(0, 'Screensize'));
end
%% 
disp('---------------------------------------------------------------------')

%% save warp if not existing
if ~exist([dataFolder,filesep,datasetName,'_flow.mat'],'file')
    v = mainSuper.v;
    warpOp = mainSuper.warpOp;
    save([dataFolder,filesep,datasetName,'_flow.mat'],'v','warpOp');
end
