% Manual Image Registration Using Control Points
% Allows the user to manually pick corresponding points between two images,
% then computes and applies a geometric transform.

clear; close all; 

%% 1) Read input and reference images
% Use specified file paths
refRGB   = imread('/Volumes/School/Thesis/data/captures/HSI/cactus_D50_pPhoto_after.png');
inputRGB   = imread('/Volumes/School/Thesis/data/captures/film_scans/icc_cactus_led_fuji_underexp_d50.png');

%% 2) Crop ROIs interactively) Crop ROIs interactively
% Crop input image
figure('Name','Crop Input','NumberTitle','off');
imshow(inputRGB); title('Draw rectangle for INPUT ROI');
hIn = drawrectangle(); wait(hIn); rectIn = round(hIn.Position);
inputRGB = imcrop(inputRGB, rectIn); close;

% Crop reference image
figure('Name','Crop Reference','NumberTitle','off');
imshow(refRGB); title('Draw rectangle for REFERENCE ROI');
hRef = drawrectangle(); wait(hRef); rectRef = round(hRef.Position);
refRGB = imcrop(refRGB, rectRef); close;

%% 3) Manual control point selection

movingPointsFile = 'movingPoints.mat';
fixedPointsFile  = 'fixedPoints.mat';

% Default to empty
movingPoints = [];
fixedPoints = [];

% Load saved points if both files exist
if exist(movingPointsFile, 'file') && exist(fixedPointsFile, 'file')
    load(movingPointsFile, 'movingPoints');
    load(fixedPointsFile,  'fixedPoints');
end

% Ensure correct format and check if both are non-empty
useInitialPoints = ~isempty(movingPoints) && ~isempty(fixedPoints) && ...
                   size(movingPoints,2) == 2 && size(fixedPoints,2) == 2;

% Launch cpselect accordingly
if useInitialPoints
    [movingPoints, fixedPoints] = cpselect(inputRGB, refRGB, ...
        movingPoints, fixedPoints, 'Wait', true);
else
    [movingPoints, fixedPoints] = cpselect(inputRGB, refRGB, 'Wait', true);
end

% Save updated points
save(movingPointsFile, 'movingPoints');
save(fixedPointsFile, 'fixedPoints');

% Cleanup: ensure only X,Y
movingPoints = double(movingPoints(:,1:2));
fixedPoints  = double(fixedPoints(:,1:2));


%% 4) Estimate geometric transform from control points from control points
% Choose transform type: 'projective', 'affine', 'similarity'
transformType = 'projective';
tform = fitgeotrans(movingPoints, fixedPoints, transformType);

%% 5) Apply transform and crop
outputView = imref2d(size(refRGB(:,:,1)));
registered = imwarp(inputRGB, tform, 'OutputView', outputView);

% Crop to overlapping region
mask = any(registered>0,3);
props = regionprops(mask,'BoundingBox');
bbox = isempty(props) * [1 1 size(refRGB,2) size(refRGB,1)] + ...
       ~isempty(props) * round(props(1).BoundingBox);
registeredCropped = imcrop(registered, bbox);
refCropped        = imcrop(refRGB,    bbox);

%% 6) Save outputs


outputDir = '/Volumes/School/Thesis/data/captures/registered';  % or specify a full path
registeredFile = fullfile(outputDir, 'cactus3_reg_fuji_led_after.png');
croppedRefFile = fullfile(outputDir, 'cactus3_ref_hsi_fuji_led_after_underexp.png');
imwrite(registeredCropped, registeredFile);
imwrite(refCropped,  croppedRefFile);

%% 7) Display results
figure('Name','Manual Registration','NumberTitle','off','Position',[100 100 1200 500]);
subplot(1,3,1), imshow(inputRGB),       title('Original Input');
subplot(1,3,2), imshow(refRGB),         title('Original Reference');
subplot(1,3,3), imshowpair(refCropped, registeredCropped), title('Manual Registered');

% Optional: show difference
figure('Name','Difference','NumberTitle','off');
imshow(imabsdiff(refCropped, registeredCropped), []);

disp('Manual registration complete. Use cpselect to refine control points if needed.');



%% SSIM & Dice
% reg = double(imread('/Volumes/School/Thesis/thesis-repo/code/matlab/img_diff/results/yoda_registered_hsi.png'));
% ref = double(imread('/Volumes/School/Thesis/thesis-repo/code/matlab/img_diff/results/yoda_reference_hsi.png'));
reg = double(registeredCropped);   % Registered input
ref = double(refCropped);   % Cropped reference

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