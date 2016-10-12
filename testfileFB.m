clear all;close all;clc;

addpath(genpath(cd));

%% load data
dataset = 'munich';
dataPath = ['..',filesep,'superResolutionShared',filesep,'data',filesep,dataset,filesep];
dataFile = [dataPath,'data.mat'];

load(dataFile);

if (~exist('imageSequenceLarge','var'))
    imageSequenceLarge = 0;
end

imageSequenceLarge = imageSequenceLarge(1:100,1:100,:);

%% init

factor = 1;
alpha = 0.01;
beta = 0.02;
numFrames = 5;


mainSuper = jointSuperResolutionFB(imageSequenceLarge(:,:,1:numFrames),alpha,beta,factor,'gtU',imageSequenceLarge);
mainSuper.verbose = 1;
%%
mainSuper.init;

mainSuper.numMainIt = 1;
%% run

mainSuper.run;

%%

[errorU,errorV] = mainSuper.calculateErrors

%%
close all;

for i=1:numFrames
    if (exist('imageSequenceLarge','var'))
        err = ssim(mainSuper.u(:,:,i),imageSequenceLarge(:,:,i));
    end
    figure(1);imagesc(mainSuper.u(:,:,i),[0,1]);drawnow; axis image;colormap(gray);pause
end

