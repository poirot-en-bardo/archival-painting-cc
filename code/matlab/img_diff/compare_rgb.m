% Three-Image SURF-Only Registration and Stitching with Subplots
% Uses SURF features, projective transform, and interactive ROI.

clear; close all; clc;

%% 1) Read input and reference images
inputRGB = imread('/Volumes/School/Thesis/data/comparison/icc_cactus_halogen_fuji_exp0_before.png');
refRGB   = imread('/Volumes/School/Thesis/data/comparison/icc_cactus_halogen_D65_interm.png');
% refRGB = imread('/Volumes/School/Thesis/data/comparison/icc_cactus_halogen_D65_pPhoto_before.png');

%% 2) Select ROIs interactively
figure('Name','Select Input ROI','NumberTitle','off');
imshow(inputRGB); title('Draw rectangle for INPUT');
h1 = drawrectangle(); rect1 = round(h1.Position);
inputCrop = imcrop(inputRGB, rect1);
close;

figure('Name','Select Reference ROI','NumberTitle','off');
imshow(refRGB); title('Draw rectangle for REFERENCE');
h2 = drawrectangle(); rect2 = round(h2.Position);
refCrop = imcrop(refRGB, rect2);
close;

%% 3) Resize crops to match
[h1,w1,~] = size(inputCrop);
[h2,w2,~] = size(refCrop);
targetSize = [min(h1,h2), min(w1,w2)];
inputRes = imresize(inputCrop, targetSize);
refRes   = imresize(refCrop,   targetSize);

%% 4) Register and crop both images
[inputReg] = register_images(inputRes, refRes);

%% 5) Show summary subplots
figure('Name','Registration Results','NumberTitle','off','Position',[100 100 1200 600]);
subplot(2,3,1), imshow(inputRes),    title('1: Input Crop');
subplot(2,3,2), imshow(refRes),      title('2: Reference Crop');
subplot(2,3,3), showMatchedFeatures(rgb2gray(inputRes), rgb2gray(refRes), ...
    inputReg.inliersMoving, inputReg.inliersFixed, 'montage');
title('3: Inlier SURF Matches');
subplot(2,3,4), imshowpair(inputReg.refCropped, inputReg.registered),
    title('4: Registered vs. Reference');
subplot(2,3,5), imshow(imabsdiff(inputReg.refCropped, inputReg.registered),[]),
    title('5: Abs Difference');


%%
outputDir = 'results';  % or specify a full path
registeredFile = fullfile(outputDir, 'registered_interm.png');
croppedRefFile = fullfile(outputDir, 'ref_interm.png');
imwrite(inputReg.registered, registeredFile);
imwrite(inputReg.refCropped,  croppedRefFile);
fprintf('Saved registered input to %s\n', registeredFile);
fprintf('Saved cropped reference to %s\n', croppedRefFile);



%% Local function for SURF registration with dual cropping
function result = register_images(inputImg, refImg)
    % Convert to grayscale and single
    grayIn  = im2single(rgb2gray(inputImg));
    grayRef = im2single(rgb2gray(refImg));

    % Detect & extract SURF features
    ptsIn  = detectSURFFeatures(grayIn,  'MetricThreshold',500);
    ptsRef = detectSURFFeatures(grayRef, 'MetricThreshold',500);
    [fIn,  vptsIn]  = extractFeatures(grayIn,  ptsIn);
    [fRef, vptsRef] = extractFeatures(grayRef, ptsRef);

    % Match features
    idxPairs = matchFeatures(fIn, fRef, 'MatchThreshold',1.5, 'MaxRatio',0.9, 'Unique',true);
    matchedMoving = vptsIn(idxPairs(:,1));
    matchedFixed  = vptsRef(idxPairs(:,2));

    % Estimate projective transform with RANSAC
    [tform, inlierIdx] = estimateGeometricTransform2D(...
        matchedMoving.Location, matchedFixed.Location, 'projective', ...
        'MaxDistance',4, 'Confidence',99.9, 'MaxNumTrials',3000);

    % Warp input image into ref coordinate frame
    outputView = imref2d(size(refImg(:,:,1)));
    regFull = imwarp(inputImg, tform, 'OutputView', outputView);

    % Find bounding box of valid (non-zero) pixels
    mask = any(regFull>0,3);
    props = regionprops(mask,'BoundingBox');
    bbox  = round(props(1).BoundingBox);

    % Crop both ref and registered to the same box
    regCropped = imcrop(regFull, bbox);
    refCropped = imcrop(refImg , bbox);

    % Prepare inlier points for montage
    inMoving = matchedMoving(inlierIdx);
    inFixed  = matchedFixed(inlierIdx);

    % Return results
    result.registered     = regCropped;
    result.refCropped     = refCropped;
    result.tform          = tform;
    result.inliersMoving  = inMoving;
    result.inliersFixed   = inFixed;
end