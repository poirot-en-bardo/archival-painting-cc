clear; close all;

% INPUTS
ref_mat = '/home/oem/eliza/mac-shared/registered/yoda/yoda_after.mat';
mov_mat = '/home/oem/eliza/mac-shared/registered/yoda/yoda_halogen_fuji_exp0_data_1.mat';
outFolder = '/home/oem/eliza/mac-shared/registered_output';
if ~exist(outFolder, 'dir'); mkdir(outFolder); end

% Load .mat files
ref = load(ref_mat);
mov = load(mov_mat);

% --- White balance the moving ProPhoto RGB image
figure; imshow(mov.RGB_img);
title('Draw rectangle around white patch and double-click');
hW = drawrectangle(); wait(hW);
whitePos = round(hW.Position);
patch_rgb = mov.RGB_img(whitePos(2):(whitePos(2)+whitePos(4)-1), whitePos(1):(whitePos(1)+whitePos(3)-1), :);
white_patch = squeeze(mean(mean(patch_rgb,1),2));
rgbWB = mov.RGB_img ./ reshape(white_patch,1,1,3);
rgbWB = max(min(rgbWB,1),0);
mov_rgb_wb = rgbWB;
%%
% Crop moving image before registration
figure; imshow(mov_rgb_wb);
title('Draw rectangle to crop moving image and double-click');
hCrop = drawrectangle(); wait(hCrop);
cropPos = round(hCrop.Position);
mov_rgb_wb = imcrop(mov_rgb_wb, cropPos);
mov.XYZ_img = imcrop(mov.XYZ_img, cropPos);

%%
% --- Resize both reference and moving images to same size
ref_img = im2uint8(ref.RGB_img);
mov_img = im2uint8(mov_rgb_wb);
% [h1,w1,~] = size(mov_img);
% [h2,w2,~] = size(ref_img);
% targetSize = [min(h1,h2), min(w1,w2)];
% mov_img = imresize(mov_img, targetSize);
% ref_img = imresize(ref_img, targetSize);
% mov.XYZ_img = imresize(mov.XYZ_img, targetSize);

% --- SURF Registration
grayMov = rgb2gray(mov_img); grayRef = rgb2gray(ref_img);
pts1 = detectSIFTFeatures(grayRef);
pts2 = detectSIFTFeatures(grayMov);
[f1, vpts1] = extractFeatures(grayRef, pts1);
[f2, vpts2] = extractFeatures(grayMov, pts2);
idxPairs = matchFeatures(f2, f1, 'MatchThreshold', 8, 'MaxRatio', 0.4, 'Unique', true);
matchedMoving = vpts2(idxPairs(:,1));
matchedFixed = vpts1(idxPairs(:,2));
[tform, inlierIdx] = estimateGeometricTransform2D(matchedMoving.Location, matchedFixed.Location, ...
    'projective', 'MaxDistance', 4, 'Confidence', 96, 'MaxNumTrials', 10000);

matchedPoints1 = matchedMoving(inlierIdx, :);  % matched1 is Nx2
matchedPoints2 = matchedFixed(inlierIdx, :);
% figure;
% showMatchedFeatures(ref_img, mov_img, matchedPoints1, matchedPoints2, 'montage');
% title('Inlier matched points');
% 
% % --- Visualize matches 
figure; showMatchedFeatures(ref_img, mov_img, matchedMoving, matchedFixed, 'montage');
title('Matched points');


% --- Warp & mask out black borders
outputView = imref2d(size(ref_img));
reg_rgb = imwarp(mov_img, tform, 'OutputView', outputView);
reg_xyz = imwarp(mov.XYZ_img, tform, 'OutputView', outputView);
reg_prophoto = imwarp(mov_rgb_wb, tform, 'OutputView', outputView);

% --- Crop to valid overlap region
mask = any(reg_rgb > 0, 3);
props = regionprops(mask, 'BoundingBox');
if isempty(props); error('No valid overlapping area found.'); end
bbox = round(props(1).BoundingBox);
reg_rgb = imcrop(reg_rgb, bbox);
reg_xyz = imcrop(reg_xyz, bbox);
reg_prophoto = imcrop(reg_prophoto, bbox);
ref_cropped = imcrop(ref.XYZ_img, bbox);
ref_rgb_prophoto = imcrop(ref.RGB_img, bbox);

% --- Visualization
figure; 
subplot(1,3,1); imshow(mov_img); title('Moving (Cropped+WB)');
subplot(1,3,2); imshow(ref_img); title('Reference');
subplot(1,3,3); imshowpair(reg_prophoto, ref_rgb_prophoto); title('Overlay (Registered)');

%%
% --- Convert to Lab and Linear RGB
reg_lab = xyz2lab(reshape(reg_xyz, [], 3));
reg_rgb_lin = max(0, min(1, xyz2rgb(reg_xyz / 100, 'ColorSpace','linear-rgb','WhitePoint','d50')));
ref_lab = xyz2lab(reshape(ref_cropped, [], 3));
ref_rgb_lin = max(0, min(1, xyz2rgb(ref_cropped / 100, 'ColorSpace','linear-rgb','WhitePoint','d50')));

[h_xyz, w_xyz, ~] = size(reg_xyz);
[h_ref, w_ref, ~] = size(ref_cropped);

reg_lab = reshape(reg_lab, h_xyz, w_xyz, 3);
ref_lab = reshape(ref_lab, h_ref, w_ref, 3);

% --- Save outputs
[~, mov_name, ~] = fileparts(mov_mat);
[~, ref_name, ~] = fileparts(ref_mat);
saveProPhotoTIFF(im2uint16(reg_prophoto), fullfile(outFolder, [mov_name '_reg_prophoto.tif']));
imwrite(im2uint16(reg_rgb_lin), fullfile(outFolder, [mov_name '_reg_linear.tif']));
saveProPhotoTIFF(im2uint16(ref_rgb_prophoto), fullfile(outFolder, [mov_name '_ref_prophoto.tif']));
imwrite(im2uint16(ref_rgb_lin), fullfile(outFolder, [ref_name '_ref_linear.tif']));

% --- Save .mat
save(fullfile(outFolder, [mov_name '_and_' ref_name '_registered.mat']), ...
    'reg_xyz', 'reg_lab', 'reg_prophoto', 'ref_cropped', 'ref_lab', 'ref_rgb_prophoto');

% --- Evaluate registration
reg_gray = mat2gray(rgb2gray(reg_prophoto));
ref_gray = mat2gray(rgb2gray(ref_rgb_prophoto));
ssimVal = ssim(reg_gray, ref_gray);
regMask = imbinarize(reg_gray);
refMask = imbinarize(ref_gray);
diceCoeff = 2 * nnz(regMask & refMask) / (nnz(regMask) + nnz(refMask));
fprintf('SSIM: %.4f, Dice: %.3f\n', ssimVal, diceCoeff);

% --- Visualization
figure; 
subplot(1,3,1); imshow(mov_img); title('Moving (Cropped+WB)');
subplot(1,3,2); imshow(ref_img); title('Reference');
subplot(1,3,3); imshowpair(reg_prophoto, ref_rgb_prophoto); title('Overlay (Registered)');
