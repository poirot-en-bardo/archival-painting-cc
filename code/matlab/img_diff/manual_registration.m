% Manual Image Registration Using Control Points
% Allows the user to manually pick corresponding points between two images,
% then computes and applies a geometric transform.

clear; close all; clc;

%% 1) Read input and reference images
% Use specified file paths
inputRGB = imread('/Volumes/School/Thesis/data/comparison/icc_cactus_halogen_D65_before.png');
refRGB   = imread('/Volumes/School/Thesis/data/comparison/icc_cactus_halogen_D65_interm.png');

%% 2) Crop ROIs interactively) Crop ROIs interactively
% Crop input image
figure('Name','Crop Input','NumberTitle','off');
imshow(inputRGB); title('Draw rectangle for INPUT ROI');
hIn = drawrectangle(); rectIn = round(hIn.Position);
inputRGB = imcrop(inputRGB, rectIn); close;

% Crop reference image
figure('Name','Crop Reference','NumberTitle','off');
imshow(refRGB); title('Draw rectangle for REFERENCE ROI');
hRef = drawrectangle(); rectRef = round(hRef.Position);
refRGB = imcrop(refRGB, rectRef); close;

%% 3) Manual control point selection
% Launch control point selection tool
movingPointsFile = 'movingPoints.mat';
fixedPointsFile  = 'fixedPoints.mat';

% If you already have saved points, load them; otherwise, run cpselect
if exist(movingPointsFile,'file') && exist(fixedPointsFile,'file')
    load(movingPointsFile,'movingPoints');
    load(fixedPointsFile,'fixedPoints');
else
    [movingPoints, fixedPoints] = cpselect(inputRGB, refRGB, 'Wait', true);
    save(movingPointsFile,'movingPoints');
    save(fixedPointsFile, 'fixedPoints');
end
% Ensure points are Nx2 numeric arrays
movingPoints = double(movingPoints);
fixedPoints = double(fixedPoints);
if size(movingPoints,2) > 2
    movingPoints = movingPoints(:,1:2);
end
if size(fixedPoints,2) > 2
    fixedPoints = fixedPoints(:,1:2);
end

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
imwrite(registeredCropped, 'results/manual_registered.png');
imwrite(refCropped,        'results/manual_reference.png');

%% 7) Display results
figure('Name','Manual Registration','NumberTitle','off','Position',[100 100 1200 500]);
subplot(1,3,1), imshow(inputRGB),       title('Original Input');
subplot(1,3,2), imshow(refRGB),         title('Original Reference');
subplot(1,3,3), imshowpair(refCropped, registeredCropped), title('Manual Registered');

% Optional: show difference
figure('Name','Difference','NumberTitle','off');
imshow(imabsdiff(refCropped, registeredCropped), []);

disp('Manual registration complete. Use cpselect to refine control points if needed.');
