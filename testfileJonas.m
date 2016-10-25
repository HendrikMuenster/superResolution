%
% Video Super Resolution
% mangled variant of Hendrik's framework, updated to jointSuperResolutionFB from 12.10
%


clearvars;
%clc;
%close all



%% Load data 
datasetName = 'city';
color = 1;
startFrame = 1;
numFrames = 5;

dataPath = ['..',filesep,'superResolutionData',filesep,datasetName,filesep];

if color
    dataFile = [dataPath,'data_color.mat'];
else
    dataFile = [dataPath,'data.mat'];
end


load(dataFile);
addpath(genpath(cd)); % only load paths if data location is correct as well



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
mainSuper = jointSuperResolutionJonas(imageSequenceStart,'gtU', imageSequenceLarge,'datasetName',datasetName);


%% Set variables

% Prodcedure
mainSuper.factor = 4;                      % magnification factor
mainSuper.numMainIt = 2;                   % number of total outer iterations
mainSuper.verbose = 1;                     % enable intermediate output
mainSuper.profiler = 0;                    % enable profiling

% Problem parameters
mainSuper.regU = 'infTGV';%'TGV';                    % regU type, 'TV' / 'TGV' (is TGV-2)
mainSuper.regUq = 1;                       % TV/TGV-2 exponent
mainSuper.regTGV = sqrt(2);                % TV-TGV weight - this value is coupled to alpha !
mainSuper.alpha = 0.005; % regU weights
mainSuper.kOpts.delta  = 0.1;              % blur l2 penalties
mainSuper.kOpts.zeta   = 0.1;              % blur Tikh penalties

mainSuper.regV = 'Huber';                  % regV type
mainSuper.beta = 0.1;                      % regV weights
mainSuper.eta  = 0.006;                     % warp weight
mainSuper.sigma = 0;%5;                       % patch size

mainSuper.gamma = 1;                       % Adaptive Coupling Term weights
mainSuper.kOpts.initKx =  exp(-(-3:3).^2 / 1.6);  % initial separable kernel (x)
mainSuper.kOpts.initKy =  exp(-(-3:3).^2 / 1.6);   % initial separable kernel (y)


%% INIT flow field and solvers
mainSuper.init;


%% Solve joint problem by alternating through u, v and k subproblems

mainSuper.run;

%% Show error margin

[errorU,errorV] = mainSuper.calculateErrors;


%% visualize video and compare to bicubic upsampling
if mainSuper.verbose
    vid = implay(mainSuper.u,2);  % persists through close all ( \_('')_/� ) - use this to compare visual quality to previous iterations
    set(vid.Parent, 'Name', [mainSuper.regU,'-',num2str(mainSuper.regUq),...
                            ' with alpha = ',num2str(mainSuper.alpha(end)),...
                            ', Huber with beta = ',num2str(mainSuper.beta(end)), ...
                            ', gamma is ',num2str(mainSuper.gamma(end))]);
    set(vid.Parent, 'Position',get(0, 'Screensize'));
    
    
    %figure(403), imagesc(imresize(mainSuper.imageSequenceSmall(:,:,3),4)),colormap('gray'),caxis([0,1]),axis image
end
%% 
disp('---------------------------------------------------------------------')
close all
%beep