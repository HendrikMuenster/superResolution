clear all;close all;clc;

addpath(genpath(cd));

%% load data
dataset = 'mapVideo';
dataPath = ['data',filesep,dataset,filesep];
dataFile = [dataPath,'map.mat'];

load(dataFile);



%% init

factor = 2;
alpha = 0.01;
beta = 0.001;
numFrames = 5;

mainSuper = jointSuperResolution(imageSequenceSmall(:,:,1:numFrames),alpha,beta,factor);
mainSuper.verbose = 1;
%%
mainSuper.init;

mainSuper.numMainIt = 1;
%% run


mainSuper.run;

%%

%[errorU,errorV] = mainSuper.calculateErrors

%%
close all;
% 
% for i=1:numFrames
%     figure(1);imagesc(mainSuper.u(:,:,i));drawnow; axis image;pause
% end

