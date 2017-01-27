%
% Video Super Resolution

%clearvars;
pPSNR = [];
pSSIM = [];
CruncherPath = matlab.desktop.editor.getActive;
cd(fileparts(CruncherPath.Filename));

%% Data properties
datasetName = 'city';
startFrame = 1;
numFrames = 13;


testCases = {'FB','F-I','I-B','FMB'};

for ii = 1:length(testCases)
    
    %% Load Video and code
    dataFolder = '/windows/DataJonas/ScieboLocalFolder/Data/videos_scenes/';
    [imageSequenceSmall,imageSequenceLarge] = LoadImSequence([dataFolder,filesep,datasetName],startFrame,numFrames);
    addpath(genpath(cd)); % only load paths if data location was error free
    
    
    
    
    
    %% Init algorithm class thing
    mainSuper = jointSuperResolutionMinimal(imageSequenceSmall);
    
    %% Set variables
    
    % Procedure
    mainSuper.factor        = 4;                   % magnification factor, remember to change
    mainSuper.numMainIt     = 1;                   % number of total outer iterations
    mainSuper.verbose       = 1;                   % enable intermediate output, 1 is text, 2 is image
    
    % Problem parameters
    mainSuper.alpha1        = 0.01;                % warp weight
    mainSuper.alpha2        = 0.01;                % TV weight
    mainSuper.beta          = 0.1;                 % regU weights
    mainSuper.kappa         = 0.5;                 % regularization pendulum value, set to NaN for standard TV, or to 1 for isoWarp-TV
    mainSuper.kOpts.delta   = 0.5;                 % blur Tikh penalties
    
    
    
    %% INIT flow fields and solve
    mainSuper.init;
    mainSuper.run;
    
    outTest = mainSuper.result1(20:end-20,20:end-20,:,:);
    groundT = imageSequenceLarge(20:end-20,20:end-20,:,:);
    psnrVal = zeros(numFrames,1); ssimVal = zeros(numFrames,1);
    for i = 1:numFrames
        psnrVal(i) = psnr(outTest(:,:,:,i),groundT(:,:,:,i));
        ssimVal(i) = ssim(outTest(:,:,:,i),groundT(:,:,:,i));
    end
    figure, imshow(outTest(:,:,:,numFrames)), title(['Last frame - ',testCases{ii}])
    figure, imshow(outTest(:,:,:,1)), title(['First frame - ',testCases{ii}])
    
    pPSNR = cat(1,pPSNR,psnrVal');
    pSSIM = cat(1,pSSIM,ssimVal');
end
%%

figure(),plot(pPSNR'),legend('FB','F-I','I-B','FMB');

figure(),plot(pSSIM'),legend('FB','F-I','I-B','FMB');
%%
disp('---------------------------------------------------------------------')
%close all
%beep