%
% Video Processor script 
% Given an input folder of ground truth images data/ground_truth,
% the script produces quantized and clipped bicubic downsampled input frames
% in folder data/input in the first step and in the second step runs the 
% MMC algorithm for all frames, writing result1, result2, u-w and w into output folders
%
% The size of the batches has to be chosen manually relative to GPU memory.

%% Control

% Choose Data
data = {'surfing','city','calendar_high','foreman'};


% Specify data location !
dataFolder = '/windows/DataJonas/ScieboLocalFolder/Data/video_scenes_long/';

% specify sampling method
downSamplingMethod = 'bicubic';
Antialiasing       = true;
factor             = 4;        % output is to be upsampled  by this factor
writeNewInput      = 0;


% MMC parameters
% 31 is ok for calendar
% 23 is ok for city 
% 89 is ok for foreman 
% 9 for surfing 
batchSizeList = [9,23,31,89]; % depends on video resolution
overlap       = 1; % batch overlap is one as one frame is fixed as described in supp. material


%% Preprocessor ( Data generator)

for ii = 1:length(data)
    
    
    % Disable write input if necessary
    if ~writeNewInput
        break
    end
    % check input data
    gt_adress = [dataFolder,data{ii},filesep,'ground_truth',filesep];
    fileStruct = [dir([gt_adress,'*.jpg']);dir([gt_adress,'*.png']); dir([gt_adress,'*.tif'])];
    numFramesTotal = length(fileStruct);
    
    
    % Create folder if it does not exist and clear previous images
    if ~exist([dataFolder,filesep,data{ii},'/input'],'dir')
        mkdir([dataFolder,filesep,data{ii},'/input'])
    else
        delete([dataFolder,filesep,data{ii},'input/*.png'])
    end
    
    
    for jj = 1:numFramesTotal
        
        imgHigh = im2double(imread([dataFolder,filesep,data{ii},filesep,'ground_truth',filesep,fileStruct(jj).name]));
        
        % unify input to always be RGB
        if size(imgHigh,3) == 1
            imgHigh = repmat(imgHigh,1,1,3);
        end
        
        % run downsampling
        imgLow = imresize(imgHigh,1/factor,downSamplingMethod,'Antialiasing',Antialiasing);
        
        % crop overshooting values
        imgLow = max(min(imgLow,1),0);
        
        % write downsampled image to folder as png, i.e. quantizing
        imwrite(imgLow,[dataFolder,filesep,data{ii},filesep,'input',filesep,'Frame',num2str(jj,'%03d'),'.png']);
        
    end
    
    % write sizes of large image
    [xLarge,yLarge,~] = size(imgHigh);
    % write sizes of small image
    [xSmall,ySmall,~] = size(imgLow);
    
    disp(['Input data written for data set ',data{ii}]);
end


%% MultiFrame Motion Coupling in connected batches

for ii = 1:length(data)
    
    % output folder name
    out_name = '/outputMMC_MON_3';
    % Read input file structure
    dataAdress = [dataFolder,filesep,data{ii},filesep,'input',filesep];
    fileStruct = dir([dataAdress,'*.png']);
    numFramesTotal = length(fileStruct);
    
    % 
    testIm  = imread([dataAdress,fileStruct(1).name]);
    [xSmall,ySmall,~] = size(testIm);
    batchSize = batchSizeList(ii);
    batchCount = ceil(numFramesTotal / (batchSize-overlap));
    %batchEnd   = numFramesTotal-(batchCount-1)*(batchSize-overlap);
    batchEnd   = batchCount*(batchSize-overlap)-numFramesTotal; 
    if batchEnd == 0 
        batchEnd = batchSize; % :>
    end
    
    % Make output folder / clear old output folder
    if ~exist([dataFolder,filesep,data{ii},out_name],'dir')
        mkdir([dataFolder,filesep,data{ii},out_name])
        mkdir([dataFolder,filesep,data{ii},out_name,'/standard'])
        mkdir([dataFolder,filesep,data{ii},out_name,'/tmpStab'])
        mkdir([dataFolder,filesep,data{ii},out_name,'/uw'])
        mkdir([dataFolder,filesep,data{ii},out_name,'/w'])
    end
    
    % Init fix frames
    u0_frame = [];
    w0_frame = [];
    for kk = 1:batchCount
        
        disp(['Running superResolution on ',data{ii},' batch ',num2str(kk),' of ',num2str(batchCount)]);
        
        if kk ~=batchCount  
            % Find ids of current batch
            fileNames = ((kk-1)*(batchSize-overlap)+1) : ((kk-1)*(batchSize-overlap)+batchSize);
            
            % Read image sequence batch
            imageSequenceSmall = zeros(xSmall,ySmall,3,batchSize);
            for jj = 1:batchSize
                imageSequenceSmall(:,:,:,jj) = im2double(imread([dataAdress,fileStruct(fileNames(jj)).name]));
            end    
        else
            % Find remaining ids
            fileNames = (numFramesTotal-batchEnd+1) : numFramesTotal; 
            
            % Read remaining image sequence
            imageSequenceSmall = zeros(xSmall,ySmall,3,batchEnd);
            for jj = 1:batchEnd
                imageSequenceSmall(:,:,:,jj) = im2double(imread([dataAdress,fileStruct(fileNames(jj)).name]));
            end
        end
        mainSuper = MultiframeMotionCoupling(imageSequenceSmall,'u0_frame',u0_frame,'w0_frame',w0_frame);
        mainSuper.factor         = factor;       % magnification factor
        mainSuper.verbose        = 1;            % enable intermediate output, 1 is text, 2 is image
        
        % Problem parameters
        mainSuper.alpha          = 0.01;                % temporal weight
        mainSuper.beta           = 0.2;                 % flow field complexity
        mainSuper.kappa          = 0.25;                % regularization pendulum
        mainSuper.flowDirection  = 'forward';

        % Operator details
        mainSuper.interpMethod   = 'average';
        mainSuper.k              = fspecial('gaussian',7,sqrt(0.6));
        % Run the thing
        mainSuper.init;
        mainSuper.run;
        
        % write output images
        for jj = 1:length(fileNames)
            out_id = fileNames;
            imwrite(mainSuper.result1(:,:,:,jj),[dataFolder,filesep,data{ii},out_name,'/standard/Frame',num2str(out_id(jj),'%03d'),'.png']);
            imwrite(mainSuper.result2(:,:,:,jj),[dataFolder,filesep,data{ii},out_name,'/tmpStab/Frame',num2str(out_id(jj),'%03d'),'.png']);
            uw = mainSuper.u(:,:,jj)-mainSuper.w(:,:,jj);
            w = mainSuper.w(:,:,jj);
            
            % align uw to first frame
            if kk== 1 && jj == 1 % first frame
                uw_min = min(uw(:)); uw_max = max(uw(:))-uw_min;
            end
            % align w to some rough bounds
            w_min = -0.25; w_max = 0.5;
            % Fit to [0,1]
            uw = (uw-uw_min)/uw_max;
            w = (w-w_min)/w_max;

            imwrite(uw,[dataFolder,filesep,data{ii},out_name,'/uw/Frame',num2str(out_id(jj),'%03d'),'.png']);
            imwrite(w,[dataFolder,filesep,data{ii},out_name,'/w/Frame',num2str(out_id(jj),'%03d'),'.png']);
        end
        
        % get fix frames for next iteration
        u0_frame = mainSuper.u(:,:,end);
        w0_frame = mainSuper.w(:,:,end);

        
        clear imageSequenceSmall
    end
end
