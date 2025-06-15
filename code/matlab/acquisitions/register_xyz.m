clear; close all;

% INPUTS
ref_mat = '/home/oem/eliza/mac-shared/registered/cactus/cactus_after.mat';
mov_mat = '/home/oem/eliza/mac-shared/registered/cactus/cactus_led_fuji_underexp.mat';
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
% --- Crop reference image safely before resizing
figure; imshow(ref.RGB_img);
title('Draw rectangle to crop reference image and double-click');
hCropRef = drawrectangle(); wait(hCropRef);
cropRefPos = round(hCropRef.Position);

% Clip to image bounds
[h_ref, w_ref, ~] = size(ref.RGB_img);
x1 = max(1, cropRefPos(1));
y1 = max(1, cropRefPos(2));
x2 = min(w_ref, x1 + cropRefPos(3) - 1);
y2 = min(h_ref, y1 + cropRefPos(4) - 1);

% Create cropped region
ref.RGB_img = ref.RGB_img(y1:y2, x1:x2, :);
ref.XYZ_img = ref.XYZ_img(y1:y2, x1:x2, :);  % Also crop XYZ if available

%%

% --- Resize both reference and moving images to same size
[h1, w1, ~] = size(ref.RGB_img);
[h2, w2, ~] = size(mov_rgb_wb);
targetSize = [min(h1,h2), min(w1,w2)];

% Resize RGB images
ref.RGB_img     = imresize(ref.RGB_img,     targetSize);
mov_rgb_wb      = imresize(mov_rgb_wb,      targetSize);

% Resize XYZ images
ref.XYZ_img     = imresize(ref.XYZ_img,     targetSize);
mov.XYZ_img     = imresize(mov.XYZ_img,     targetSize);


%%
ref_img = im2uint8(ref.RGB_img);
mov_img = im2uint8(mov_rgb_wb);
mov_img = im2uint8(mov_rgb_wb * 1.7);  % Increase brightness by 50%
mov_img = min(mov_img, 255);           % Clip to valid uint8 range
% mov_img = im2uint8(imadjust(mov_rgb_wb, [], [], 0.7));  % Gamma < 1 brightens image

% --- SURF Registration
grayMov = rgb2gray(mov_img); 
grayRef = rgb2gray(ref_img);
% pts1 = detectSURFFeatures(grayRef, 'MetricThreshold', 1000);
% pts2 = detectSURFFeatures(grayMov, 'MetricThreshold', 1000);
pts1 = detectSIFTFeatures(grayRef);
pts2 = detectSIFTFeatures(grayMov);
[f1, vpts1] = extractFeatures(grayRef, pts1);
[f2, vpts2] = extractFeatures(grayMov, pts2);
idxPairs = matchFeatures(f2, f1, 'MatchThreshold', 6, 'MaxRatio', 0.3, 'Unique', true);
matchedMoving = vpts2(idxPairs(:,1));
matchedFixed = vpts1(idxPairs(:,2));
[tform, inlierIdx] = estimateGeometricTransform2D(matchedMoving.Location, matchedFixed.Location, ...
    'projective', 'MaxDistance', 7, 'Confidence', 95, 'MaxNumTrials', 5000);

matchedPoints1 = matchedMoving(inlierIdx, :);  % matched1 is Nx2
matchedPoints2 = matchedFixed(inlierIdx, :);
% figure;
% showMatchedFeatures(ref_img, mov_img, matchedPoints1, matchedPoints2, 'montage');
% title('Inlier matched points');

% 
% --- Visualize matches 
figure; showMatchedFeatures(ref_img, mov_img, matchedMoving, matchedFixed, 'montage');
title('Matched points');
%
% manual

% mov_img = im2uint8(mov_rgb_wb * 1.1);  % Increase brightness by 50%
% mov_img = min(mov_img, 255); 
% 
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


%
% --- Warp into reference coordinate frame
outputView = imref2d(size(ref_img));
reg_rgb       = imwarp(mov_img,     tform, 'OutputView', outputView);
reg_xyz       = imwarp(mov.XYZ_img, tform, 'OutputView', outputView);
reg_prophoto  = imwarp(mov_rgb_wb,  tform, 'OutputView', outputView);

% --- Create a valid mask from the registered RGB
mask = any(reg_rgb > 0, 3);
props = regionprops(mask, 'BoundingBox');
if isempty(props)
    error('No valid overlapping area found after warping.');
end

% --- Crop all images using the same bounding box
bbox = round(props(1).BoundingBox);
reg_rgb         = imcrop(reg_rgb, bbox);
reg_xyz         = imcrop(reg_xyz, bbox);
reg_prophoto    = imcrop(reg_prophoto, bbox);
ref_prophoto = imcrop(ref.RGB_img, bbox);
ref_cropped     = imcrop(ref.XYZ_img, bbox);

% --- Final size check
[h1, w1, ~] = size(reg_prophoto);
[h2, w2, ~] = size(ref_prophoto);
if h1 ~= h2 || w1 ~= w2
    minH = min(h1, h2);
    minW = min(w1, w2);
    reg_prophoto = imresize(reg_prophoto, [minH, minW]);
    ref_prophoto = imresize(ref_prophoto, [minH, minW]);
    reg_rgb = imresize(reg_rgb, [minH, minW]);
    reg_xyz = imresize(reg_xyz, [minH, minW]);
    ref_cropped = imresize(ref_cropped, [minH, minW]);
end

% --- Visualization
figure; 
subplot(2,2,1); imshow(mov_img); title('Moving (Cropped+WB)');
subplot(2,2,2); imshow(ref_img); title('Reference');
subplot(2,2,3); imshowpair(reg_prophoto, ref_prophoto); title('Overlay (Registered)');
subplot(2,2,4); imshow(imabsdiff(im2uint8(reg_prophoto), im2uint8(ref_prophoto)), []);
title('Absolute Difference');

%
% --- Convert to Lab and Linear RGB
reg_lab = xyz2lab(reshape(reg_xyz, [], 3));
reg_rgb_lin = max(0, min(1, xyz2rgb(reg_xyz / 100, 'ColorSpace','linear-rgb','WhitePoint','d50')));
ref_lab = xyz2lab(reshape(ref_cropped, [], 3));
ref_rgb_lin = max(0, min(1, xyz2rgb(ref_cropped / 100, 'ColorSpace','linear-rgb','WhitePoint','d50')));

[h_xyz, w_xyz, ~] = size(reg_xyz);
[h_ref, w_ref, ~] = size(ref_cropped);

reg_lab = reshape(reg_lab, h_xyz, w_xyz, 3);
ref_lab = reshape(ref_lab, h_ref, w_ref, 3);


% --- Evaluate registration
reg_gray = mat2gray(rgb2gray(reg_prophoto));
ref_gray = mat2gray(rgb2gray(ref_prophoto));
ssimVal = ssim(reg_gray, ref_gray);
regMask = imbinarize(reg_gray);
refMask = imbinarize(ref_gray);
diceCoeff = 2 * nnz(regMask & refMask) / (nnz(regMask) + nnz(refMask));
fprintf('SSIM: %.4f, Dice: %.3f\n', ssimVal, diceCoeff);
%%
% --- Save outputs
[~, mov_name, ~] = fileparts(mov_mat);
[~, ref_name, ~] = fileparts(ref_mat);

% --- Create subfolder for this moving image
mov_out_dir = fullfile(outFolder, mov_name);
if ~exist(mov_out_dir, 'dir'); mkdir(mov_out_dir); end

% --- Save outputs in subfolder
saveProPhotoTIFF(im2uint16(reg_prophoto), fullfile(mov_out_dir, [mov_name '_reg_prophoto.tif']));
imwrite(im2uint16(reg_rgb_lin), fullfile(mov_out_dir, [mov_name '_reg_linear.tif']));
saveProPhotoTIFF(im2uint16(ref_prophoto), fullfile(mov_out_dir, [mov_name '_ref_prophoto.tif']));
imwrite(im2uint16(ref_rgb_lin), fullfile(mov_out_dir, [ref_name '_ref_linear.tif']));

% --- Save .mat file in subfolder
save(fullfile(mov_out_dir, [mov_name '_and_' ref_name '_registered.mat']), ...
    'reg_xyz', 'reg_lab', 'reg_prophoto', 'ref_cropped', 'ref_lab', 'ref_prophoto');
%%
% --- Visualization
figure; 
subplot(1,3,1); imshow(mov_rgb_wb); title('Moving (Film)');
subplot(1,3,2); imshow(ref_img); title('Reference (Scan)');
%%
subplot(1,3,3); imshowpair(reg_prophoto, ref_prophoto); title('Overlay Registered');


%%