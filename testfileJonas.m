%
% Video Super Resolution
% mangled variant of Hendrik's framework, updated to jointSuperResolutionFB from 12.10
%


clearvars;
%clc;
close all


%% Load data 
datasetName = 'city';
color = 0;

dataPath = ['..',filesep,'superResolutionData',filesep,datasetName,filesep];

if color
    dataFile = [dataPath,'data_color.mat']; %#ok<UNRCH>
else
    dataFile = [dataPath,'data.mat'];
end

try
    load(dataFile);
catch 
    cd ..\
    load(dataFile);
end

numFrames = 5;

%% Init algorithm class thing
mainSuper = jointSuperResolutionJonas(imageSequenceSmall(:,:,1:numFrames),'gtU',imageSequenceLarge,'datasetName',datasetName);


%% Set variables

% Prodcedure
mainSuper.factor = 4;                      % magnification factor
mainSuper.numMainIt = 1;                   % number of total outer iterations
mainSuper.verbose = 1;                     % enable intermediate output
mainSuper.profiler = 0;                    % enable profiling

% Problem parameters
mainSuper.regU = 'TGV';                    % regU type, 'TV' / 'TGV' (is TGV-2)
mainSuper.regUq = 0.5;                     % TV/TGV-2 exponent
mainSuper.regTGV = 1;                    % TV-TGV weight
mainSuper.alpha = 0.001;%[0.1,0.01,0.001]; % regU weights
mainSuper.kOpts.delta  = 0.05;              % blur l2 penalties
mainSuper.kOpts.zeta   = 0.1;             % blur Tikh penalties

mainSuper.regV = 'Huber';                  % regV type
mainSuper.beta = 0.05;%[0.3,0.2,0.1]; % regV weights

mainSuper.gamma = 1;                       % Coupling Term weights



%% INIT flow field and solvers
mainSuper.init;


%% Solve joint problem by alternating through u, v and k subproblems

mainSuper.run;

%% Show error margin

[errorU,errorV] = mainSuper.calculateErrors;

%% visualize video and compare to bicubic upsampling
if mainSuper.verbose
    implay(mainSuper.u,2);  % persists through close all ( \_('')_/¯ ) - use this to compare visual quality to previous iterations
    figure(403), imagesc(imresize(mainSuper.imageSequenceSmall(:,:,3),4)),colormap('gray'),caxis([0,1]),axis image
end
%% 
disp('---------------------------------------------------------------------')
