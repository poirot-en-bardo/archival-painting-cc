% Three-Image SURF-Only Registration and Stitching with Subplots
% Uses SURF features, projective transform, and interactive ROI.

clear; close all; clc;

%% 1) Read input and reference images
refRGB   = imread('/Volumes/School/Thesis/data/captures/HSI/cactus_D50_pPhoto_after.png');
inputRGB   = imread('/Volumes/School/Thesis/data/captures/film_scans/icc_cactus_led_fuji_underexp_d50.png');

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

%%
% INPUTS
ref_mat = '/home/oem/eliza/mac-shared/registered/yoda/yoda_after.mat';
mov_mat = '/home/oem/eliza/mac-shared/registered/yoda/yoda_halogen_fuji_exp0_data_1.mat';
outFolder = '/home/oem/eliza/mac-shared/registered_output';
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% Load .mat files
ref = load(ref_mat);
mov = load(mov_mat);

% --- White balance the moving ProPhoto RGB image
figure; imshow(mov.RGB_img);
title('Draw rectangle around white patch and double-click');
hW = drawrectangle();
wait(hW);
whitePos = hW.Position;
xw = max(1, floor(whitePos(1)));
yw = max(1, floor(whitePos(2)));
ww = floor(whitePos(3));
hw = floor(whitePos(4));
xw_end = min(size(mov.RGB_img,2), xw + ww - 1);
yw_end = min(size(mov.RGB_img,1), yw + hw - 1);
patch_rgb = mov.RGB_img(yw:yw_end, xw:xw_end, :);
white_patch = squeeze(mean(mean(patch_rgb,1),2));

% Apply white balance
rgbWB = mov.RGB_img ./ reshape(white_patch,1,1,3);
rgbWB = max(min(rgbWB,1),0);
mov_rgb_wb = rgbWB;
%%
% Crop moving image before registration
figure; imshow(mov_rgb_wb);
title('Draw rectangle to crop moving image and double-click');
hCrop = drawrectangle();
wait(hCrop);
cropPos = round(hCrop.Position);
xc = max(1, cropPos(1));
yc = max(1, cropPos(2));
wc = cropPos(3);
hc = cropPos(4);
x_end = min(size(mov_rgb_wb, 2), xc + wc - 1);
y_end = min(size(mov_rgb_wb, 1), yc + hc - 1);

mov_rgb_wb = mov_rgb_wb(yc:y_end, xc:x_end, :);
mov.XYZ_img = mov.XYZ_img(yc:y_end, xc:x_end, :);


% Use white-balanced RGB for registration
ref_img = im2uint8(ref.RGB_img);
mov_img = im2uint8(mov_rgb_wb);
%%
inputCrop = mov_img;
refCrop = ref_img;
%% 3) Resize crops to match
[h1,w1,~] = size(inputCrop);
[h2,w2,~] = size(refCrop);
targetSize = [min(h1,h2), min(w1,w2)];
inputRes = imresize(inputCrop, targetSize);
refRes   = imresize(refCrop,   targetSize);

%% Register and crop both images
[inputReg] = register_images(inputRes, refRes);

%% Show summary subplots
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
outputDir = '/Volumes/School/Thesis/data/captures/registered';  % or specify a full path
registeredFile = fullfile(outputDir, 'cactus3_reg_fuji_led_after.png');
croppedRefFile = fullfile(outputDir, 'cactus3_ref_hsi_fuji_led_after_underexp.png');
imwrite(inputReg.registered, registeredFile);
imwrite(inputReg.refCropped,  croppedRefFile);
fprintf('Saved registered input to %s\n', registeredFile);
fprintf('Saved cropped reference to %s\n', croppedRefFile);

%% Evaluating the registration: SSIM & Dice

reg = double(inputReg.registered);   % Registered input
ref = double(inputReg.refCropped);   % Cropped reference

% Convert to double in [0,1] before binarizing
reg_norm = mat2gray(rgb2gray(double(reg)));
ref_norm = mat2gray(rgb2gray(double(ref)));

regMask = imbinarize(reg_norm);
refMask = imbinarize(ref_norm);

dice_coeff = 2 * nnz(regMask & refMask) / (nnz(regMask) + nnz(refMask));
fprintf('Dice Overlap Coefficient = %.3f\n', dice_coeff);

reg_gray = mat2gray(rgb2gray(reg));  % Normalize to [0,1]
ref_gray = mat2gray(rgb2gray(ref));
ssimVal = ssim(reg_gray, ref_gray);
fprintf('SSIM (grayscale) = %.4f\n', ssimVal);

%% Local function for SIFT registration with dual cropping
function result = register_images(inputImg, refImg)
    % Convert to grayscale and single
    grayIn  = im2single(rgb2gray(inputImg));
    grayRef = im2single(rgb2gray(refImg));

    ptsIn  = detectSIFTFeatures(grayIn);
    ptsRef = detectSIFTFeatures(grayRef);
    [fIn,  vptsIn]  = extractFeatures(grayIn,  ptsIn);
    [fRef, vptsRef] = extractFeatures(grayRef, ptsRef);

    % Match features with Lowe's ratio test
    idxPairs = matchFeatures(fIn, fRef, ...
        'MatchThreshold', 5, ...
        'MaxRatio', 0.8, ...
        'Unique', true);

    matchedMoving = vptsIn(idxPairs(:,1));
    matchedFixed  = vptsRef(idxPairs(:,2));

    % Estimate projective transform using RANSAC
    [tform, inlierIdx] = estimateGeometricTransform2D(...
        matchedMoving.Location, matchedFixed.Location, ...
        'projective', 'MaxDistance', 7, 'Confidence', 95, 'MaxNumTrials', 5000);

    % Warp input image into reference coordinate frame
    outputView = imref2d(size(refImg(:,:,1)));
    regFull = imwarp(inputImg, tform, 'OutputView', outputView);

    % Mask out black borders from warping
    mask = any(regFull > 0, 3);
    props = regionprops(mask, 'BoundingBox');
    bbox  = round(props(1).BoundingBox);

    % Crop both images to shared area
    regCropped = imcrop(regFull, bbox);
    refCropped = imcrop(refImg, bbox);

    % Return results
    result.registered     = regCropped;
    result.refCropped     = refCropped;
    result.tform          = tform;
    result.inliersMoving  = matchedMoving(inlierIdx);
    result.inliersFixed   = matchedFixed(inlierIdx);
end


%% manual registration


% movingPointsFile = 'movingPoints.mat';
% fixedPointsFile = 'fixedPoints.mat';
% movingPoints = []; fixedPoints = [];
% 
% if exist(movingPointsFile, 'file') && exist(fixedPointsFile, 'file')
%     load(movingPointsFile, 'movingPoints');
%     load(fixedPointsFile,  'fixedPoints');
% end
% 
% useInitialPoints = ~isempty(movingPoints) && ~isempty(fixedPoints) && ...
%                    size(movingPoints,2) == 2 && size(fixedPoints,2) == 2;
% 
% if useInitialPoints
%     [movingPoints, fixedPoints] = cpselect(mov_img, ref_img, movingPoints, fixedPoints, 'Wait', true);
% else
%     [movingPoints, fixedPoints] = cpselect(mov_img, ref_img, 'Wait', true);
% end
% 
% save(movingPointsFile, 'movingPoints');
% save(fixedPointsFile,  'fixedPoints');
% 
% movingPoints = double(movingPoints(:,1:2));
% fixedPoints  = double(fixedPoints(:,1:2));
% 
% % --- Compute transform and warp
% tform = fitgeotrans(movingPoints, fixedPoints, 'projective');