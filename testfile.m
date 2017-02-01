%
% Video Super Resolution

addpath(genpath(cd));

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
%dataFolder = 'C:\Users\Hendrik\Dropbox\Uni\Projects\2016 - SuperResolutionMunich\superResolutionData\';
[imageSequenceSmall,imageSequenceLarge] = LoadImSequence([dataFolder,filesep,datasetName],startFrame,numFrames,factor);    



%% Load precomputed flow field if existing
saveAndLoad = 0;
if exist([dataFolder,filesep,datasetName,'_flow.mat'],'file') && saveAndLoad
    load([dataFolder,filesep,datasetName,'_flow.mat']);
    figure(2), imagesc(flowToColorV2(squeeze(v(:,:,cslice,:)))); axis image; drawnow
else
    v = 0;
    warpOp = 0;
end



%% Init algorithm class thing
mainSuper = MultiframeMotionCoupling(imageSequenceSmall,'flowField',v,'warpOp',warpOp);

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
mainSuper.interpAA     = false;              % The antialiasing switch is complicated in 2017a
kernel = @(x) double(abs(x)<=2)/4;
width  = 7;
mainSuper.interpKernel = {kernel,width};

% "Motion" blur
mainSuper.k = fspecial('gaussian',7,sqrt(0.6)); 

            

%% INIT flow field and solvers
mainSuper.init;


%% Solve super resolution problem in u
tic
mainSuper.run;
toc

%% Show error margin of central patch

outImage = mainSuper.result1(20:end-20,20:end-20,:,ceil(numFrames/2));
psnrErr = psnr(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2)));
ssimErr = ssim(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2)));
disp(['PSNR (central patch, central slice): ',num2str(psnrErr),' dB']);
disp(['SSIM (central patch, central slice): ',num2str(ssimErr),' ']);

%% visualize video
if mainSuper.verbose > 0
    vid = implay(mainSuper.result1,2);  % persists through close all ( \_('')_/ ) - use this to compare visual quality to previous iterations
    set(vid.Parent, 'Name', ['u of InfAddTV with kappa = ',num2str(mainSuper.kappa),...
                            ' and alphas ',num2str(mainSuper.alpha1),', ',num2str(mainSuper.alpha2),...
                            ', regV is Huber with beta = ',num2str(mainSuper.beta), ...
                            ' - all for ',num2str(mainSuper.numMainIt),' iterations']);
    set(vid.Parent, 'Position',get(0, 'Screensize'));
    
%     vid = implay(mainSuper.result2,2);  % persists through close all ( \_('')_/ ) - use this to compare visual quality to previous iterations
%     set(vid.Parent, 'Name', ['u-w of InfAddTV with kappa = ',num2str(mainSuper.kappa),...
%                             ' and alphas ',num2str(mainSuper.alpha1),', ',num2str(mainSuper.alpha2),...
%                             ', regV is Huber with beta = ',num2str(mainSuper.beta), ...
%                             ' - all for ',num2str(mainSuper.numMainIt),' iterations']);
%     set(vid.Parent, 'Position',get(0, 'Screensize'));
end
%%
% %% write central image to file
% fileNaming = ['results',filesep,datasetName,'_TVinfAdd_u -  alpha', ...
%     num2str(mainSuper.alpha1,4),', kappa',num2str(mainSuper.kappa,4),', its ',num2str(mainSuper.numMainIt),'.png'];
% imwrite(mainSuper.result1(:,:,:,cslice),fileNaming);
% %% write central image to file
% fileNaming = ['results',filesep,datasetName,'_TVinfAdd_u-w -  alpha1', ...
%     num2str(mainSuper.alpha1,4),', alpha2',num2str(mainSuper.alpha2,4),', its ',num2str(mainSuper.numMainIt),'.png'];
% imwrite(mainSuper.result2(:,:,:,cslice),fileNaming);

%% 
disp('---------------------------------------------------------------------')

%% save warp if not existing
if ~exist([dataFolder,filesep,datasetName,'_flow.mat'],'file') && saveAndLoad 
    v = mainSuper.v;
    warpOp = mainSuper.warpOp;
    save([dataFolder,filesep,datasetName,'_flow.mat'],'v','warpOp');
end
