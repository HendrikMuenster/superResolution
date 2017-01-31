%
% Video Super Resolution
close all
clearvars;
pPSNR = [];
pSSIM = [];
CruncherPath = matlab.desktop.editor.getActive;
cd(fileparts(CruncherPath.Filename));

%% Data properties
datasetName = 'walk_high';
startFrame = 1;
numFrames = 13;


testCases = {'FB','F-I','I-B','FMB'};

for ii = 1:length(testCases)
    
    %% Load Video and code
    dataFolder = '/windows/DataJonas/ScieboLocalFolder/Data/videos_scenes/';
    [imageSequenceSmall,imageSequenceLarge] = LoadImSequence([dataFolder,filesep,datasetName],startFrame,numFrames);
    addpath(genpath(cd)); % only load paths if data location was error free
    
    
    
    
    
    %% Init algorithm class thing
    mainSuper = jointSuperResolutionMinimal(imageSequenceSmall,'testCase',testCases{ii});
    
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
    mainSuper.interpAA     = false;
    width = 7;
    kernel = @(x) double(abs(x)<=2)/factor; % this is the old average kernel
    mainSuper.interpKernel = {kernel,width};

    % "Motion" blur:
    mainSuper.k = fspecial('gaussian',7,sqrt(0.6));
    
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

figure(),plot(pPSNR','-x'),legend('FB','F-I','I-B','FMB');title('PSNR values per frame');

figure(),plot(pSSIM','-x'),legend('FB','F-I','I-B','FMB');title('SSIM values per frame');

drawnow
%%
disp('---------------------------------------------------------------------')
%close all
%beep