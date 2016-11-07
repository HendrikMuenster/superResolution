%
% Video Super Resolution
% mangled variant of Hendrik's framework, updated to jointSuperResolutionFB from 12.10
%


%clearvars;
%clc;
%close all



%% Load data 
datasetName = 'city';
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

% Add noise
% for ii = 1:numFrames
%     imageSequenceSmall(:,:,:,ii) = imnoise(imageSequenceSmall(:,:,:,ii), 'salt & pepper',0.1); %#ok<SAGROW>
% end
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
%mainSuper = jointSuperResolutionJonas(imageSequenceStart,'gtU', imageSequenceLarge,'datasetName',datasetName);
mainSuper = jointSuperResolutionJonas(imageSequenceStart,'gtU',imageSequenceLarge);

%% Set variables

% Prodcedure
mainSuper.factor = 4;                      % magnification factor
mainSuper.numMainIt = 1;                   % number of total outer iterations
mainSuper.verbose = 1;                     % enable intermediate output, 1 is text, 2 is image
mainSuper.profiler = 0;                    % enable profiling

% Problem parameters
mainSuper.regU = 'infAdd';                  % regU type, 'TV' / 'TGV' (is TGV-2)
mainSuper.regUq = 1;                       % TV/TGV-2 exponent
mainSuper.regTGV = sqrt(2);                % TV-TGV weight - this value is coupled to alpha !
mainSuper.alpha = 0.01;                  % regU weights
mainSuper.kOpts.delta  = 0;                % blur l2 penalties
mainSuper.kOpts.zeta   = 0.01;              % blur Tikh penalties

mainSuper.regV = 'Huber';                  % regV type
mainSuper.beta = 0.1;                     % regV weights
mainSuper.eta  = 0.01;                    % warp weight
mainSuper.opts.nsize = 0;                  % radius of local boundary, choose 0 for no local boundaries
mainSuper.opts.offset = 0.1;               % offset of local boundary
mainSuper.gamma = 1/eps;                   % outlier removal, to to 1/eps for total outlier removal, but 1 is usually enough
mainSuper.sigma = 0;                       % patchwise extension of outlier removal                   

mainSuper.kOpts.initKx =  exp(-(-3:3).^2 / 1.2);  % initial separable kernel (x)
mainSuper.kOpts.initKy =  exp(-(-3:3).^2 / 1.2);  % initial separable kernel (y)


%% INIT flow field and solvers
mainSuper.init;


%% Solve joint problem by alternating through u, v and k subproblems

mainSuper.run;


%% Show error margin

[errorU,errorV] = mainSuper.calculateErrors;


%% visualize video and compare to bicubic upsampling
if mainSuper.verbose > 0
    vid = implay(mainSuper.u,2);  % persists through close all ( \_('')_/¯ ) - use this to compare visual quality to previous iterations
    set(vid.Parent, 'Name', [mainSuper.regU,'-',num2str(mainSuper.regUq),...
                            ' with alpha = ',num2str(mainSuper.alpha(end)),...
                            ', Huber with beta = ',num2str(mainSuper.beta(end)), ...
                            ', eta is ',num2str(mainSuper.eta(end))]);
    set(vid.Parent, 'Position',get(0, 'Screensize'));
    
    
    %figure(403), imagesc(imresize(mainSuper.imageSequenceSmall(:,:,3),4)),colormap('gray'),caxis([0,1]),axis image
end

%% write central image to file
fileNaming = ['results',filesep,datasetName,'_',mainSuper.regU,'-',num2str(mainSuper.regUq),' alpha', ...
    num2str(mainSuper.alpha(end),4),' eta',num2str(mainSuper.eta,4),' its ',num2str(mainSuper.numMainIt),'.png'];
cslice = floor(numFrames/2);
imwrite(mainSuper.u(:,:,:,cslice),fileNaming);

%% 
disp('---------------------------------------------------------------------')
%close all
%beep