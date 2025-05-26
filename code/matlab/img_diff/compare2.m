% Three-Image Combined Registration and Advanced Refinement
% Uses SURF + KAZE detection, numeric point matching, projective transform, interactive ROI,
% and multi-stage intensity-based refinement (affine then translation).

clear; close all; clc;

%% 1) Read input and reference images
% inputRGB = imread('/Volumes/School/Thesis/data/comparison/icc_cactus_halogen_fuji_exp0_before.png');
% refRGB   = imread('/Volumes/School/Thesis/data/comparison/icc_cactus_halogen_D65_interm.png');
inputRGB = imread('registered_input_interm.png');
inputRGB = imread('cropped_reference_interm.png');
refRGB = imread('registered_input_interm.png');


%% 2) Select ROIs interactively
figure('Name','Select Input ROI','NumberTitle','off');
imshow(inputRGB); title('Draw rectangle for INPUT');
h1 = drawrectangle(); rect1 = round(h1.Position);
inputCrop = imcrop(inputRGB, rect1); close;

figure('Name','Select Reference ROI','NumberTitle','off');
imshow(refRGB); title('Draw rectangle for REFERENCE');
h2 = drawrectangle(); rect2 = round(h2.Position);
refCrop = imcrop(refRGB, rect2); close;

%% 3) Resize crops to match
[h1,w1,~] = size(inputCrop);
[h2,w2,~] = size(refCrop);
targetSize = [min(h1,h2), min(w1,w2)];
inputRes = imresize(inputCrop, targetSize);
refRes   = imresize(refCrop,   targetSize);

%% 4) Register and crop both images
[inputReg] = register_and_advanced_refine(inputRes, refRes);

%% 4.1) Save outputs for Fourier difference later
outputDir = pwd;
imwrite(inputReg.registered, fullfile(outputDir, 'registered_input_interm.png'));
imwrite(inputReg.refCropped, fullfile(outputDir, 'cropped_reference_interm.png'));

%% 5) Show summary subplots
figure('Name','Registration Results','NumberTitle','off','Position',[100 100 1200 600]);
subplot(2,3,1), imshow(inputRes),    title('1: Input Crop');
subplot(2,3,2), imshow(refRes),      title('2: Reference Crop');
subplot(2,3,3), showMatchedFeatures(rgb2gray(inputRes), rgb2gray(refRes), ...
    inputReg.inMoving, inputReg.inFixed, 'montage');
title('3: Inlier Matches');
subplot(2,3,4), imshowpair(inputReg.refCropped, inputReg.registered),
    title('4: Registered vs. Reference');
subplot(2,3,5), imshow(imabsdiff(inputReg.refCropped, inputReg.registered),[]),
    title('5: Abs Difference');

%% Local function: Combined detectors + multi-stage intensity-based refinement
function result = register_and_advanced_refine(inputImg, refImg)
    % Convert to grayscale
    grayIn  = im2single(rgb2gray(inputImg));
    grayRef = im2single(rgb2gray(refImg));

    % 1) Feature detection (SURF + KAZE)
    ptsSurfIn  = detectSURFFeatures(grayIn,  'MetricThreshold',50).selectStrongest(5000);
    ptsSurfRef = detectSURFFeatures(grayRef, 'MetricThreshold',50).selectStrongest(5000);
    ptsKazeIn  = detectKAZEFeatures(grayIn).selectStrongest(5000);
    ptsKazeRef = detectKAZEFeatures(grayRef).selectStrongest(5000);

    % Extract features
    [fSurfIn,  vSurfIn]  = extractFeatures(grayIn,  ptsSurfIn,  'Upright', true);
    [fSurfRef, vSurfRef] = extractFeatures(grayRef, ptsSurfRef, 'Upright', true);
    [fKazeIn,  vKazeIn]  = extractFeatures(grayIn,  ptsKazeIn,  'Upright', true);
    [fKazeRef, vKazeRef] = extractFeatures(grayRef, ptsKazeRef, 'Upright', true);

    % Match features numerically
    idxSurf = matchFeatures(fSurfIn, fSurfRef, 'MaxRatio',0.6, 'Unique',true, 'MatchThreshold', 40);
    surfMovLoc = ptsSurfIn.Location(idxSurf(:,1),:);
    surfFixLoc = ptsSurfRef.Location(idxSurf(:,2),:);
    idxKaze = matchFeatures(fKazeIn, fKazeRef, 'MaxRatio',0.6, 'Unique',true, 'MatchThreshold', 40);
    kazeMovLoc = ptsKazeIn.Location(idxKaze(:,1),:);
    kazeFixLoc = ptsKazeRef.Location(idxKaze(:,2),:);

    % Combine both sets
    allMov = [surfMovLoc; kazeMovLoc];
    allFix = [surfFixLoc; kazeFixLoc];

    % 2) Initial projective RANSAC
    [tformProj, inlierIdx] = estimateGeometricTransform2D(...
        allMov, allFix, 'projective', 'MaxDistance',3, 'Confidence',99.99, 'MaxNumTrials',20000);

    % Warp and crop
    outputView = imref2d(size(refImg(:,:,1)));
    reg1 = imwarp(inputImg, tformProj, 'OutputView', outputView);
    mask = any(reg1>0,3);
    props = regionprops(mask,'BoundingBox');
    if isempty(props)
        bbox = [1,1,size(refImg,2),size(refImg,1)];
    else
        bbox = round(props(1).BoundingBox);
    end
    regCropped = imcrop(reg1, bbox);
    refCropped = imcrop(refImg,  bbox);

    % 3) Affine intensity-based refinement
    [optA, metA] = imregconfig('multimodal');
    tformAff = imregtform(rgb2gray(regCropped), rgb2gray(refCropped), 'affine', optA, metA);
    regAff = imwarp(regCropped, tformAff, 'OutputView', imref2d(size(refCropped)));

    % 4) Final translation refinement
    [optT, metT] = imregconfig('multimodal');
    tformTrans = imregtform(rgb2gray(regAff), rgb2gray(refCropped), 'translation', optT, metT);
    regFinal = imwarp(regAff, tformTrans, 'OutputView', imref2d(size(refCropped)));

    % Gather inliers for display (use original RANSAC inliers)
    inMoving = allMov(inlierIdx, :);
    inFixed  = allFix(inlierIdx, :);

    % Return
    result.registered     = regFinal;
    result.refCropped     = refCropped;
    result.inMoving       = inMoving;
    result.inFixed        = inFixed;
    result.tformProj      = tformProj;
    result.tformAff       = tformAff;
    result.tformTrans     = tformTrans;
end
