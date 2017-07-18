%% MULTIFRAME MOTION COUPLING FOR VIDEO SUPER RESOLUTION
%
%
% Be sure to initialize the submodules and to compile them
% following their instructions


clearvars;
%addpath(genpath(cd)); use floated flexBox



%% Data properties
data = {'tube3','city','foreman','walk','foliage','calendar','surfer','surferdog','penguin','temple','sheets','wave'};
dataFolder = '/windows/DataJonas/ScieboLocalFolder/Data/videos_scenes/';
writeFolder = '/windows/DataJonas/ScieboLocalFolder/Data/MMC_result/';
startFrame = 1;
numFramesList = [ones(6,1)*13;ones(6,1)*5];
factor  = 4;             % Magnification factor
vl_setupnn();
%% Run the thing
for kk = 1:length(data)
    disp_('Running on dataset',data{kk}, '.........')
    numFrames = numFramesList(kk);
    
    % Create folder if it does not exist and clear previous images
    if ~exist([writeFolder,filesep,data{kk},'/input'],'dir')
        mkdir([writeFolder,filesep,data{kk},'/input'])
    else
        delete([writeFolder,filesep,data{kk},'input/*.png'])
    end
    
    % Load images
    [imageSequenceSmall,imageSequenceLarge] = LoadImSequence([dataFolder,filesep,data{kk}],startFrame,numFrames,factor,'bicubic');
    
    
    
    %% Construct algorithm object
    % Input: RGB-Time matlab array
    t1 = tic;
    mainSuper = MultiframeMotionCoupling(imageSequenceSmall);
    
    %% Set variables (these are the standard parameters)
    
    % Procedure
    mainSuper.factor        = factor;              % magnification factor
    mainSuper.verbose       = 1;                   % enable intermediate output, 1 is text, 2 is image
    mainSuper.framework     = 'flexBox';           % Choose framework for super resolution problem
    
    % Problem parameters
    mainSuper.alpha         = 0.01;                % regularizer weight
    mainSuper.beta          = 0.2;                 % flow field complexity
    mainSuper.kappa         = 0.25;                % regularization pendulum
    mainSuper.flowDirection = 'forward';           % flow field direction
    
    % Operator details
    mainSuper.interpMethod  = 'average';            % Downsampling operator D
    mainSuper.k = fspecial('gaussian',7,sqrt(0.6));% Blur operator B
    mainSuper.VDSRFlag      = true;                 % Enable VDSR initial guess
    
    
    %% Init flow field and solvers and run MMC
    t2 = tic;
    mainSuper.init_v0;
    OFTime(kk) = toc(t2); %#ok<*SAGROW>
    
    mainSuper.init;
    t3 = tic;
    mainSuper.run;
    
    algTime(kk) = toc(t3);
    totalTime(kk) = toc(t1);
    
    %% Central point error
    
    outImage = mainSuper.result1(20:end-20,20:end-20,:,ceil(numFrames/2));
    psnrErrMid(kk) = round(psnr(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2))),2);
    ssimErrMid(kk) = round(ssim(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2))),3);
    disp(['PSNR (central patch, central slice): ',num2str(psnrErrMid(kk)),' dB']);
    disp(['SSIM (central patch, central slice): ',num2str(ssimErrMid(kk)),' ']);
    %% Average error
    psnrErrFrames = zeros(numFrames,1);
    ssimErrFrames = zeros(numFrames,1);
    for ii = 1:numFrames
        outImage = mainSuper.result1(20:end-20,20:end-20,:,ii);
        psnrErrFrames(ii) = psnr(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ii));
        ssimErrFrames(ii) = ssim(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ii));
    end
    psnrErrMean(kk) = mean(psnrErrFrames);
    ssimErrMean(kk) = mean(ssimErrFrames);
    
    
    %% Write output
    % Create folder if it does not exist and clear previous images
    if ~exist([writeFolder,filesep,data{kk},'/output'],'dir')
        mkdir([writeFolder,filesep,data{kk},'/output'])
    else
        delete([writeFolder,filesep,data{kk},'output/*.png'])
    end
    for ii = 1:numFrames
        imwrite(mainSuper.result1(:,:,:,ii),[writeFolder,filesep,data{ii},'/output/','Frame',num2str(jj,'%03d'),'.png']);
    end
    
end
