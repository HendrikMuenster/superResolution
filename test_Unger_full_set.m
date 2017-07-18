%
% Video Super Resolution

clearvars;


%% Data properties
dataset = {'tube3','city','foreman','walk','foliage','calendar','surfer','surferdog','penguin','temple','sheets','wave'};

startFrame = 1;
numFramesList = [13*ones(6,1),5*ones(6,1)];

factor = 4;                         % magnification factor


for kk = 1:length(dataset)
    
    datasetName = dataset{kk};
    numFrames = numFramesList(kk);
    cslice = ceil(numFrames/2);
    disp_('Running on',datasetName,'...')
    %% Load Video and code
    dataFolder = '/windows/DataJonas/ScieboLocalFolder/Data/videos_scenes/';
    [imageSequenceSmall,imageSequenceLarge] = LoadImSequence([dataFolder,filesep,datasetName],startFrame,numFrames,factor,'bicubic',0);
    
    
    
    %% Run Single frame motion coupling
    alpha = 0.05;
    beta = 0.1;
    
    imgSR = singleframeMotionSR(imageSequenceSmall,factor,alpha,beta);
    
    %% Show error margin
    
    outImage = imgSR(20:end-20,20:end-20,:);
    psnrErr(kk) = psnr(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2))); %#ok<*SAGROW>
    ssimErr(kk) = ssim(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ceil(numFrames/2)));
    disp(['PSNR : ',num2str(psnrErr(kk)),' dB']);
    disp(['SSIM : ',num2str(ssimErr(kk)),' ']);
    
    figure(1), imshow(outImage); title(['PSNR: ', num2str(psnrErr)]); axis image
    drawnow;
    disp('---------------------------------------------------------------------')
    
end


