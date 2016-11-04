%
% Video Super Resolution

clearvars;


%% Load data 
datasetName = 'tubeZoom2';
color = 1;
startFrame = 5;
numFrames = 5;

dataPath = ['..',filesep,'superResolutionData',filesep,datasetName,filesep];

if color
    dataFile = [dataPath,'data_color.mat'];
else
    dataFile = [dataPath,'data.mat']; %#ok<*UNRCH>
end


load(dataFile);
addpath(genpath(cd)); % only load paths if data location was error free



if color
    imageSequenceStart = im2double(imageSequenceSmall(:,:,:,startFrame:startFrame+numFrames-1));
    try 
        imageSequenceLarge = im2double(imageSequenceLarge(:,:,:,startFrame:startFrame+numFrames-1));
    catch
        imageSequenceLarge = 0;
    end
else 
    imageSequenceStart = im2double(imageSequenceSmall(:,:,startFrame:startFrame+numFrames-1));
    try
         imageSequenceLarge = im2double(imageSequenceLarge(:,:,startFrame:startFrame+numFrames-1));
    catch 
        imageSequenceLarge = 0;
    end
end



%% Init algorithm class thing
mainSuper = jointSuperResolutionMinimal(imageSequenceStart,'gtU',imageSequenceLarge);

%% Set variables

% Prodcedure
mainSuper.factor        = 4;                   % magnification factor
mainSuper.numMainIt     = 2;                   % number of total outer iterations
mainSuper.verbose       = 1;                   % enable intermediate output, 1 is text, 2 is image

% Problem parameters
mainSuper.alpha1        = 0.01;                % regU weights
mainSuper.alpha2        = 0.01;                % regU weights
mainSuper.beta          = 0.2;                 % regU weights
mainSuper.kappa         = 0.1;                 % regularization pendulum value
mainSuper.kOpts.delta   = 0.01;                % blur Tikh penalties
                  


%% INIT flow field and solvers
mainSuper.init;


%% Solve joint problem by alternating through u, v and k subproblems

mainSuper.run;


%% Show error margin

[psnrErr,ssimErr, psnrV] = mainSuper.calculateErrors;


%% visualize video
if mainSuper.verbose > 0
    vid = implay(mainSuper.result2,2);  % persists through close all ( \_('')_/¯ ) - use this to compare visual quality to previous iterations
    set(vid.Parent, 'Name', [mainSuper.regU,'-',num2str(mainSuper.regUq),...
                            ' with alpha = ',num2str(mainSuper.alpha(end)),...
                            ', Huber with beta = ',num2str(mainSuper.beta(end)), ...
                            ', eta is ',num2str(mainSuper.eta(end))]);
    set(vid.Parent, 'Position',get(0, 'Screensize'));
    
end

%% write central image to file
fileNaming = ['results',filesep,datasetName,'_',mainSuper.regU,'-',num2str(mainSuper.regUq),' alpha', ...
    num2str(mainSuper.alpha(end),4),' eta',num2str(mainSuper.eta,4),' its ',num2str(mainSuper.numMainIt),'.png'];
cslice = floor(numFrames/2);
imwrite(mainSuper.result2(:,:,:,cslice),fileNaming);

%% 
disp('---------------------------------------------------------------------')
%close all
%beep