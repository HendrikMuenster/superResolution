%% MULTIFRAME MOTION COUPLING FOR VIDEO SUPER RESOLUTION
%
%
% Be sure to initialize the submodules and to compile them
% following their instructions


clearvars;
%addpath(genpath(cd)); use floated flexBox
if exist('statblock_sr.mat','file')
    load('statblock_sr.mat');
    startHere = statblock.iterDone+1;
else
    startHere = 1;
end


%% Data properties
data = {'tube3','city','calendar_high','foliage_high','walk_high','foreman','temple','penguins','sheets','surfer','wave','surferdog'};
dataFolder = '/windows/DataJonas/ScieboLocalFolder/Data/videos_scenes/';
writeFolder = '/windows/DataJonas/ScieboLocalFolder/Data/MMC_result';
startFrame = 1;
numFramesList = [ones(6,1)*13;ones(6,1)*5];
factor  = 4;             % Magnification factor

%% Run the thing
for kk = startHere:length(data)
    disp_('Running on dataset',data{kk}, '.........')
    numFrames = numFramesList(kk);
    
    %     % Create folder if it does not exist and clear previous images
    %     if ~exist([writeFolder,filesep,data{kk},'/input'],'dir')
    %         mkdir([writeFolder,filesep,data{kk},'/input'])
    %     else
    %         delete([writeFolder,filesep,data{kk},'input/*.png'])
    %     end
    
    % Load images
    [imageSequenceSmall,imageSequenceLarge] = LoadImSequence([dataFolder,filesep,data{kk}],startFrame,numFrames,factor,'bicubic');
    %     % Write input
    %     for ii = 1:numFrames
    %         imwrite(imageSequenceSmall(:,:,:,ii),[writeFolder,filesep,data{kk},'/input/','Frame',num2str(ii,'%03d'),'.png']);
    %     end
    
    
    %% Run SingleFrame SuperResolution
    t1 = tic;
    alpha = 0.05;
    beta = 0.2;
    imgSR = singleframeMotionSR_mitzel(imageSequenceSmall,factor,alpha,beta);
    %imgSR = singleframeMotionSR_unger(imageSequenceSmall,factor,alpha,beta);
    
    totalTime(kk) = toc(t1); %#ok<*SAGROW>
    
    %% Central point error
    
    outImage = imgSR(20:end-20,20:end-20,:);
    psnrErrMid(kk) = round(psnr(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2))),2);
    ssimErrMid(kk) = round(ssim(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2))),3);
    disp(['PSNR (central patch, central slice): ',num2str(psnrErrMid(kk)),' dB']);
    disp(['SSIM (central patch, central slice): ',num2str(ssimErrMid(kk)),' ']);
    disp(totalTime)
    %% Average error
    %     psnrErrFrames = zeros(numFrames,1);
    %     ssimErrFrames = zeros(numFrames,1);
    %     for ii = 1:numFrames
    %         outImage = mainSuper.result1(20:end-20,20:end-20,:,ii);
    %         psnrErrFrames(ii) = psnr(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ii));
    %         ssimErrFrames(ii) = ssim(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ii));
    %     end
    %     psnrErrMean(kk) = mean(psnrErrFrames);
    %     ssimErrMean(kk) = mean(ssimErrFrames);
    
    
    %% Write output
    % Create folder if it does not exist and clear previous images
    %     if ~exist([writeFolder,filesep,data{kk},'/output'],'dir')
    %         mkdir([writeFolder,filesep,data{kk},'/output'])
    %     else
    %         delete([writeFolder,filesep,data{kk},'output/*.png'])
    %     end
    %     for ii = 1:numFrames
    %         imwrite(mainSuper.result1(:,:,:,ii),[writeFolder,filesep,data{kk},'/output/','Frame',num2str(ii,'%03d'),'.png']);
    %     end
    
    %% Stat Structure
    statblock.timeMidFrame(kk) = totalTime(kk);
    
    statblock.psnrErrMid(kk) = psnrErrMid(kk);
    statblock.ssimErrMid(kk) = ssimErrMid(kk);
    statblock.iterDone = kk;
    save('statblock_sr.mat','statblock');
end
