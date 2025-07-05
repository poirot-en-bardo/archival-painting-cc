clear; close all;

% ---- File paths ----
after_mat   = '/home/oem/eliza/mac-shared/registered_output/yoda_halogen_fuji_exp0/yoda_halogen_fuji_exp0_and_yoda_after_registered.mat';    % contains ref_cropped, ref_prophoto, reg_xyz, reg_prophoto
before_mat  = '/home/oem/eliza/mac-shared/registered_output/yoda_halogen_fuji_exp0_before/yoda_halogen_fuji_exp0_and_yoda_before_registered.mat';            % contains reg_xyz, reg_prophoto (before correction)
output_mat  = '/home/oem/eliza/mac-shared/registered_output/yoda_halogen_fuji_exp0_all3/yoda_halogen_fuji_exp0_all3.mat';      % output file

% ---- Load data ----
after_data  = load(after_mat);   % 'ref_cropped', 'ref_prophoto', 'reg_xyz', 'reg_prophoto'
before_data = load(before_mat);  % 'reg_xyz', 'reg_prophoto'

% ---- Assign variables for clarity ----
after_xyz   = after_data.ref_cropped;     % reference scan, XYZ
after_rgb   = after_data.ref_prophoto;    % reference scan, RGB (ProPhoto)
film_xyz    = after_data.reg_xyz;         % registered film, XYZ
film_rgb    = after_data.reg_prophoto;    % registered film, RGB (ProPhoto)

before_xyz  = before_data.ref_cropped;        % before image, XYZ
before_rgb  = before_data.ref_prophoto;   % before image, RGB (ProPhoto)
%%
% ---- Register "before" to "after" using feature-based registration ----
% Convert RGB to uint8 for SURF (can adjust method if needed)
ref_img   = im2uint8(after_rgb);
mov_img   = im2uint8(before_rgb);

% Convert to grayscale for feature extraction
gray_ref  = rgb2gray(ref_img);
gray_mov  = rgb2gray(mov_img);

% Detect and match SURF features
% pts_ref   = detectSURFFeatures(gray_ref, 'MetricThreshold', 1000);
% pts_mov   = detectSURFFeatures(gray_mov, 'MetricThreshold', 1000);
pts_ref   = detectSIFTFeatures(gray_ref);
pts_mov   = detectSIFTFeatures(gray_mov);

[feat_ref, vpts_ref] = extractFeatures(gray_ref, pts_ref);
[feat_mov, vpts_mov] = extractFeatures(gray_mov, pts_mov);

indexPairs = matchFeatures(feat_mov, feat_ref, 'Unique', true, 'MaxRatio', 0.6, 'MatchThreshold', 60);

if size(indexPairs,1) < 8
    error('Not enough matches found for reliable registration.');
end

matched_mov = vpts_mov(indexPairs(:,1));
matched_ref = vpts_ref(indexPairs(:,2));

figure('Name','Matched Feature Points');
showMatchedFeatures(ref_img, mov_img, matched_ref, matched_mov, 'montage');
title('Inlier Matched Feature Points (After vs Before)', 'FontWeight', 'bold');


[tform, inlierIdx] = estimateGeometricTransform2D(matched_mov.Location, matched_ref.Location, ...
    'similarity', 'MaxDistance', 4, Confidence=99);

%% manual 


movingPointsFile = 'movingPoints.mat';
fixedPointsFile = 'fixedPoints.mat';
movingPoints = []; fixedPoints = [];

if exist(movingPointsFile, 'file') && exist(fixedPointsFile, 'file')
    load(movingPointsFile, 'movingPoints');
    load(fixedPointsFile,  'fixedPoints');
end

useInitialPoints = ~isempty(movingPoints) && ~isempty(fixedPoints) && ...
                   size(movingPoints,2) == 2 && size(fixedPoints,2) == 2;

if useInitialPoints
    [movingPoints, fixedPoints] = cpselect(mov_img, ref_img, movingPoints, fixedPoints, 'Wait', true);
else
    [movingPoints, fixedPoints] = cpselect(mov_img, ref_img, 'Wait', true);
end

save(movingPointsFile, 'movingPoints');
save(fixedPointsFile,  'fixedPoints');

movingPoints = double(movingPoints(:,1:2));
fixedPoints  = double(fixedPoints(:,1:2));

% --- Compute transform and warp
tform = fitgeotrans(movingPoints, fixedPoints, 'projective');

% ---- Warp "before" image (XYZ and RGB) ----
outputView = imref2d(size(ref_img));
before_xyz_reg = imwarp(before_xyz, tform, 'OutputView', outputView, 'Interp', 'bilinear');
before_rgb_reg = imwarp(before_rgb, tform, 'OutputView', outputView, 'Interp', 'bilinear');
%%
figure('Name','Overlay after Registration');
imshowpair(before_rgb_reg, after_rgb); % Or use after_rgb_crop if cropped
title('Registered "Before" (overlayed on After)', 'FontWeight', 'bold');
%
figure('Name','Absolute Difference (Registered)');
imshow(imabsdiff(im2uint8(before_rgb_reg), im2uint8(after_rgb)), []);
title('Absolute Difference: Registered "Before" vs "After"', 'FontWeight', 'bold');
%%


% ---- Crop all three images to the largest common valid area ----
% Valid masks: nonzero pixels in all 3 channels
mask_after   = all(after_xyz > 0, 3);
mask_film    = all(film_xyz  > 0, 3);
mask_before  = all(before_xyz_reg > 0, 3);

mask_all     = mask_after & mask_film & mask_before;

% Find bounding box of common region
props = regionprops(mask_all, 'BoundingBox');
if isempty(props)
    error('No common valid area across all three images!');
end
bbox = round(props(1).BoundingBox);
x1 = max(1, bbox(1));
y1 = max(1, bbox(2));
x2 = min(size(after_xyz,2), x1 + bbox(3) - 1);
y2 = min(size(after_xyz,1), y1 + bbox(4) - 1);

% Crop all channels
after_xyz     = after_xyz(y1:y2, x1:x2, :);
film_xyz      = film_xyz(y1:y2, x1:x2, :);
before_xyz    = before_xyz_reg(y1:y2, x1:x2, :);

after_rgb     = after_rgb(y1:y2, x1:x2, :);
film_rgb      = film_rgb(y1:y2, x1:x2, :);
before_rgb    = before_rgb_reg(y1:y2, x1:x2, :);

% ---- Compute Lab for all images (XYZ must be HxWx3, scale if needed) ----
[h, w, ~] = size(after_xyz);
after_lab   = reshape(xyz2lab(reshape(after_xyz,[],3)), h,w,3);
film_lab    = reshape(xyz2lab(reshape(film_xyz,[],3)), h,w,3);
before_lab  = reshape(xyz2lab(reshape(before_xyz,[],3)), h,w,3);

% ---- Save to MAT file ----
save(output_mat, ...
    'after_xyz',  'after_rgb',  'after_lab', ...
    'film_xyz',   'film_rgb',   'film_lab', ...
    'before_xyz', 'before_rgb', 'before_lab', ...
    '-v7.3');

fprintf('Saved perfectly aligned images to: %s\n', output_mat);

%% ----  Visual check ----
figure('Name','Visual Check: After/Film/Before');
tiledlayout(1,3);
nexttile; imshow(before_rgb);  title('Before', 'FontWeight','bold', 'FontSize', 18);
nexttile; imshow(film_rgb);    title('Film', 'FontWeight','bold', 'FontSize', 18);
nexttile; imshow(after_rgb);   title('After ', 'FontWeight','bold', 'FontSize', 18);

