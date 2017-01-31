function [imageSequenceSmall,imageSequenceLarge] = LoadImSequence(datastr,startFrame,numFrames,factor,downSamplingMethod,beReasonable)
% load  image sequence from folder datastr
% downsample, trim and return

if nargin < 6
    beReasonable = 'ok';
end

if nargin < 5
    downSamplingMethod = 'bicubic';
end
if nargin < 4
    factor = 4;
end

olddir = cd;
cd(datastr);
fileStruct = [dir('*.jpg');dir('*.png');dir('*.tif')];

%test image
imtest = im2double(imread(fileStruct(startFrame).name));
[xLarge,yLarge,color] = size(imtest);

if strcmp(beReasonable,'ok')
    factest = 0.5;
    while xLarge > 1280 || yLarge > 1280
        xLarge = xLarge*factest;
        yLarge = yLarge*factest;
        factest = factest*0.5;
    end
end

if mod(xLarge/factor,1) || mod(yLarge/factor,1)
    error('input image is not easily dividable by target factor, do something...');
end
    
imageSequenceLarge = zeros(xLarge,yLarge,3,numFrames);
for jj = startFrame:numFrames+startFrame-1
    imageSequenceLarge(:,:,:,jj-startFrame+1) = imresize(im2double(imread(fileStruct(jj).name)),[xLarge,yLarge]);
end

% pad greyscale data
if size(imageSequenceLarge,3) == 1
    imageSequenceLarge = permute(repmat(permute(imageSequenceLarge,[1,2,4,3]),1,1,1,3),[1,2,4,3]);
end


% create downsampled data
imageSequenceSmall = zeros(xLarge/factor,yLarge/factor,color,numFrames);
for jj = 1:numFrames
    imageSequenceSmall(:,:,:,jj) = imresize(imageSequenceLarge(:,:,:,jj),1/factor,downSamplingMethod); 
end

% Correct out-of-bounds values after sampling
imageSequenceSmall(imageSequenceSmall>1) = 1;
imageSequenceSmall(imageSequenceSmall<0) = 0;

% return to old
cd(olddir);


end

