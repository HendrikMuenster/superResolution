function [imageSequenceInput] = LoadImSequenceReal(datastr,startFrame,numFrames)
% load  image sequence from folder datastr


fileStruct = [dir([datastr,filesep,'*.png']);dir([datastr,filesep,'*.jpg']);dir([datastr,filesep,'*.tif'])];

if isempty(fileStruct)
    error('invalid location path, data not found');
end

if nargin < 3
    numFrames = length(fileStruct);
end
if nargin < 2
    startFrame = 1;
end

%test image
imtest = im2double(imread([datastr,filesep,fileStruct(startFrame).name]));
[xLarge,yLarge,numColors] = size(imtest);


    
imageSequenceInput = zeros(xLarge,yLarge,numColors,numFrames);
for jj = startFrame:numFrames+startFrame-1
    if numColors == 3
        imageSequenceInput(:,:,:,jj-startFrame+1)  = im2double(imread([datastr,filesep,fileStruct(jj).name]));
    else
        imageSequenceInput(:,:,jj-startFrame+1)  = im2double(imread([datastr,filesep,fileStruct(jj).name]));
    end
end


end

