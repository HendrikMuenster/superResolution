
% compare OF fields
figure,imagesc((squeeze(v(:,:,11,1)))),caxis([0,25])
figure,imagesc((squeeze(mainSuper.v(:,:,11,1)))),caxis([0,25])

figure,imagesc((squeeze(v(:,:,11,2)))),caxis([0,25])
figure,imagesc((squeeze(mainSuper.v(:,:,11,2)))),caxis([0,25])
%% compute YLarge
for j = 1:numFrames
    imTemp = rgb2ycbcr(imageSequenceLarge(:,:,:,j));
    YLarge(:,:,j) = imTemp(:,:,1); %#ok<*SAGROW>
end

% warp YLarge
warped = mainSuper.warpOp*YLarge(:);
warpVid = reshape(warped,size(YLarge));

% results should be zero, but ...
for i = 1:numFrames
    figure,imagesc(warpVid(:,:,i)),
end

%% warp imageSeqLarge

for j = 1:3
    vecIm = squeeze(imageSequenceLarge(:,:,j,:));
    warpj = mainSuper.warpOp*vecIm(:);
    tmpRe = permute(reshape(warpj,size(vecIm)),[1,2,4,3]);
    warpVidRGB(:,:,j,:) = tmpRe;
end

%% results should be zero, but ...
for i = numFrames
    figure,imagesc(warpVidRGB(:,:,:,i)),
end


%% more testing
[X,Y] = meshgrid(1:nx*factor,1:ny*factor);
im1 = imageSequenceLarge(:,:,:,1);
im2 = imageSequenceLarge(:,:,:,2);
for j = 1:3
    warpIm(:,:,j) = interp2(im2(:,:,j),X+v(:,:,1,1),Y+v(:,:,2,1));
end
figure(5), imagesc(sum(abs(warpIm-im1),3));

%% more more testing

singleField = squeeze(v(:,:,1,:));
Wi   = warpingOperator([ny*4,nx*4],singleField);
% Remove out-of range warps
%marker = sum(abs(Wi),2) == 0;
%Wi(marker > 0,:) = 0;

im1 = imageSequenceLarge(:,:,:,1);
im2 = imageSequenceLarge(:,:,:,2);
for j = 1:3
    imj = im2(:,:,j);
    warpIm(:,:,j) = reshape(Wi*imj(:),size(imj));
end
figure(5), imagesc(sum(abs(warpIm-im1),3));

%% even more

im1 = imageSequenceLarge(:,:,:,1);
im2 = imageSequenceLarge(:,:,:,2);
for j = 1:3
    warpIm(:,:,j) = imwarp(im2(:,:,j),squeeze(v(:,:,1,:)),'cubic');
end
figure(5), imagesc(sum(abs(warpIm-im1),3)),axis image

%these are actually occlusions in the dataset, see training/occlusions/bandage_1 !!!!!