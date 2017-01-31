%
% Video Super Resolution

clearvars;
CruncherPath = matlab.desktop.editor.getActive;
cd(fileparts(CruncherPath.Filename));

%% Data properties
datasetName = 'city';
startFrame = 1;
numFrames = 13;
cslice = ceil(numFrames/2);
factor = 4;                 % magnification factor

%% Load Video and code
dataFolder = '/windows/DataJonas/ScieboLocalFolder/Data/videos_scenes/';
[imageSequenceSmall,imageSequenceLarge] = LoadImSequence([dataFolder,filesep,datasetName],startFrame,numFrames,factor);    
addpath(genpath(cd)); % only load paths if data location was error free


%% Load precomputed flow field if existing
if exist([dataFolder,filesep,datasetName,'_flow.mat'],'file') && 1
    load([dataFolder,filesep,datasetName,'_flow.mat']);
    figure(2), imagesc(flowToColorV2(squeeze(v(:,:,cslice,:)))); axis image; drawnow
else
    v = 0;
    warpOp = 0;
end



%% Init algorithm class thing
mainSuper = jointSuperResolutionMinimal(imageSequenceSmall,'flowField',v,'warpOp',warpOp);

%% Set variables

% Procedure
mainSuper.factor        = factor;              % magnification factor
mainSuper.numMainIt     = 1;                   % number of total outer iterations
mainSuper.verbose       = 1;                   % enable intermediate output, 1 is text, 2 is image

% Problem parameters
mainSuper.alpha1        = 0.01;                % temporal weight
mainSuper.alpha2        = 0.01;                % spatial weight
mainSuper.beta          = 0.1;                 % flow field complexity
mainSuper.kappa         = 0.5;                 % regularization pendulum

  
% Downsampling details
mainSuper.interpMethod = 'custom';
mainSuper.interpAA     = false;                % Increases the scaling of the kernel to prevent aliasing - MATLAB inverse, handle with care(!)

if strcmp(mainSuper.interpMethod,'custom')
    width = 7;
    kernel = @(x) double(abs(x)<=2)/factor; % this is the old average kernel
    mainSuper.interpKernel = {kernel,width};
end

% "Motion" blur:
mainSuper.k = fspecial('gaussian',7,sqrt(0.6));

%% Plot kernel comparisons
h = 0.0001;
intervalThing = -width:h:width;
sigSTD = sqrt(1.2/2);
gaussThing = @(x) 1./sqrt(2*pi*sigSTD^2).*exp(-x.^2 / 2/sigSTD^2);
figure(1), plot(intervalThing,conv(gaussThing(intervalThing),kernel(intervalThing),'same')*h);
hold on, plot(intervalThing,gaussThing(intervalThing))
plot(intervalThing,kernel(intervalThing)),

sigSTD = 1.6;
gaussThing = @(x) 1./sqrt(2*pi*sigSTD^2).*exp(-x.^2 / 2/sigSTD^2);
plot(intervalThing,gaussThing(intervalThing));
legend('Our kernel','Gaussian sigma=0.776 ','Block sampling','Gaussian sigma=1.6 ')


%% INIT flow field and solvers
mainSuper.init;


%% Solve joint problem by alternating through u, v and k subproblems

mainSuper.run;


%% Show error margin

%[psnrErr,ssimErr, psnrV] = mainSuper.calculateErrors;
outImage = mainSuper.result1(20:end-20,20:end-20,:,ceil(numFrames/2));
psnrErr = psnr(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2)));
ssimErr = ssim(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2)));
disp(['PSNR (central patch, central slice): ',num2str(psnrErr),' dB']);
disp(['SSIM (central patch, central slice): ',num2str(ssimErr),' ']);

%% visualize video
if mainSuper.verbose > 0
    vid = implay(mainSuper.result1,2);  % persists through close all ( \_('')_/� ) - use this to compare visual quality to previous iterations
    set(vid.Parent, 'Name', ['u of InfAddTV with kappa = ',num2str(mainSuper.kappa),...
                            ' and alphas ',num2str(mainSuper.alpha1),', ',num2str(mainSuper.alpha2),...
                            ', regV is Huber with beta = ',num2str(mainSuper.beta), ...
                            ' - all for ',num2str(mainSuper.numMainIt),' iterations']);
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
