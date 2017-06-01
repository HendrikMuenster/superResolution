function [imageSequenceSmall,imageSequenceLarge] = LoadImSequenceForward(datastr,startFrame,numFrames,factor,interpMethod,kernel)
% load  image sequence from folder datastr
% downsample, trim and return


olddir = cd;
cd(datastr);
fileStruct = [dir('*.jpg');dir('*.png');dir('*.tif')];

%test image
imtest = im2double(imread(fileStruct(startFrame).name));
[mLarge,nLarge,color] = size(imtest);

if mod(mLarge/factor,1) || mod(nLarge/factor,1)
    error('input image is not easily dividable by target factor, do something...');
end

% Process ground truth data to double
imageSequenceLarge = zeros(mLarge,nLarge,3,numFrames);
for jj = startFrame:numFrames+startFrame-1
    imageSequenceLarge(:,:,:,jj-startFrame+1) = im2double(imread(fileStruct(jj).name));
end

% pad greyscale data
if size(imageSequenceLarge,3) == 1
    imageSequenceLarge = permute(repmat(permute(imageSequenceLarge,[1,2,4,3]),1,1,1,3),[1,2,4,3]);
end


% create downsampled data
imageSequenceSmall = zeros(mLarge/factor,nLarge/factor,color,numFrames);

% operators:
dsOp = samplingOperator([mLarge,nLarge],[mLarge/factor,nLarge/factor],interpMethod,[],false);
dsOp = dsOp*RepConvMtx(kernel,[mLarge,nLarge]);

for jj = 1:numFrames
    imTemp = rgb2ycbcr(imageSequenceLarge(:,:,:,jj)); 
    
    imY    = imTemp(:,:,1);
    imCb_s = imresize(imTemp(:,:,2),1/factor);
    imCr_s = imresize(imTemp(:,:,3),1/factor);
    imY_s  = reshape(dsOp*imY(:),[mLarge/factor,nLarge/factor]);
                     
    imageSequenceSmall(:,:,:,jj) = ycbcr2rgb(cat(3,imY_s,imCb_s,imCr_s));                     
end

% return to old
cd(olddir);


end

