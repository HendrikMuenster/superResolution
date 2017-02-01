%
% Video Super Resolution

addpath(genpath(cd));

clearvars;
CruncherPath = matlab.desktop.editor.getActive;
cd(fileparts(CruncherPath.Filename));

%% Data properties
datasetName = 'city';
startFrame = 1;
numFrames = 3;
cslice = ceil(numFrames/2);
factor = 4;                 % magnification factor

%% Load Video and code
dataFolder = '/windows/DataJonas/ScieboLocalFolder/Data/videos_scenes/';
dataFolder = 'C:\Users\Hendrik\Dropbox\Uni\Projects\2016 - SuperResolutionMunich\superResolutionData\';
[imageSequenceSmall,imageSequenceLarge] = LoadImSequence([dataFolder,filesep,datasetName],startFrame,numFrames,factor);    
addpath(genpath(cd)); % only load paths if data location was error free


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
mainSuper.interpMethod = 'bilinear';
mainSuper.interpAA     = false;              % The antialiasing switch increases the scaling of the kernel to prevent aliasing

% Custom kernel, if necessary:
% if strcmp(mainSuper.interpMethod,'custom')
%     kernel = @(x) (sin(pi*x) .* sin(pi*x/4) + eps) ./ ((pi^2 * x.^2 / 4) + eps).*exp(-x.^2/(2*gsig^2))./sqrt(2*pi*gsig^2); % Lancz/Gauss Kernel
%     %kernel = @(x) double(abs(x)<=0.5); % this is bilinear with no aa
%     
%     width = 7;
%     figure(1), plot(-width:0.01:width,kernel(-width:0.01:width)); title('interpolation kernel'), drawnow
%     
%     mainSuper.interpKernel = {kernel,width};
% end

% "Motion" blur
mainSuper.k = fspecial('gaussian',7,1.2^2); 

            

%% INIT flow field and solvers
mainSuper.init;


%% Solve super resolution problem in u

mainSuper.run;


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
