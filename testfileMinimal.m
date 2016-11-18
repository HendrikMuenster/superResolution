%
% Video Super Resolution

clearvars;
CruncherPath = matlab.desktop.editor.getActive;
cd(fileparts(CruncherPath.Filename));

%% Data properties
datasetName = 'city';
startFrame = 1;
numFrames = 5;

%% Load Video and code
dataFolder = '/windows/DataJonas/ScieboLocalFolder/Data/videos_scenes/';
[imageSequenceSmall,imageSequenceLarge] = LoadImSequence([dataFolder,filesep,datasetName],startFrame,numFrames);    
addpath(genpath(cd)); % only load paths if data location was error free





%% Init algorithm class thing
mainSuper = jointSuperResolutionMinimal(imageSequenceSmall,'gtU',imageSequenceLarge);

%% Set variables

% Procedure
mainSuper.factor        = 4;                   % magnification factor, remember to change
mainSuper.numMainIt     = 1;                   % number of total outer iterations
mainSuper.verbose       = 1;                   % enable intermediate output, 1 is text, 2 is image

% Problem parameters
mainSuper.alpha1        = 0.01;                % regU weights
mainSuper.alpha2        = 0.01;                % regU weights
mainSuper.beta          = 0.1;                 % regU weights
mainSuper.kappa         = 0.5;                % regularization pendulum value
mainSuper.kOpts.delta   = 0.5;                 % blur Tikh penalties
                  


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
    
    vid = implay(mainSuper.result2,2);  % persists through close all ( \_('')_/� ) - use this to compare visual quality to previous iterations
    set(vid.Parent, 'Name', ['u-w of InfAddTV with kappa = ',num2str(mainSuper.kappa),...
                            ' and alphas ',num2str(mainSuper.alpha1),', ',num2str(mainSuper.alpha2),...
                            ', regV is Huber with beta = ',num2str(mainSuper.beta), ...
                            ' - all for ',num2str(mainSuper.numMainIt),' iterations']);
    set(vid.Parent, 'Position',get(0, 'Screensize'));
end
%%
%% write central image to file
fileNaming = ['results',filesep,datasetName,'_TVinfAdd_u -  alpha', ...
    num2str(mainSuper.alpha1,4),', kappa',num2str(mainSuper.kappa,4),', its ',num2str(mainSuper.numMainIt),'.png'];
cslice = ceil(numFrames/2);
imwrite(mainSuper.result1(:,:,:,cslice),fileNaming);
%% write central image to file
fileNaming = ['results',filesep,datasetName,'_TVinfAdd_u-w -  alpha1', ...
    num2str(mainSuper.alpha1,4),', alpha2',num2str(mainSuper.alpha2,4),', its ',num2str(mainSuper.numMainIt),'.png'];
cslice = ceil(numFrames/2);
imwrite(mainSuper.result2(:,:,:,cslice),fileNaming);

%% 
disp('---------------------------------------------------------------------')
%close all
%beep