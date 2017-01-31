function [imageSequenceInput] = LoadImSequenceReal(datastr,startFrame,numFrames)
% load  image sequence from folder datastr


fileStruct = dir([datastr,filesep,'*.png']);

if isempty(fileStruct)
    error('invalid location path, data not found');
end


%test image
imtest = im2double(imread([datastr,filesep,fileStruct(startFrame).name]));
[xLarge,yLarge,~] = size(imtest);


    
imageSequenceInput = zeros(xLarge,yLarge,3,numFrames);
for jj = startFrame:numFrames+startFrame-1
    imageSequenceInput(:,:,:,jj-startFrame+1)  = im2double(imread([datastr,filesep,fileStruct(jj).name]));
end


end

