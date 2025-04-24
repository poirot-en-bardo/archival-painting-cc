hdr1 = "/Volumes/School/Thesis/Data/HySpex/cactusHalogen_001_VNIR_1800_SN00841_HSNR1_24998us_2025-04-15T111238_raw.hdr";
hyspex1 = "/Volumes/School/Thesis/Data/HySpex/cactusHalogen_001_VNIR_1800_SN00841_HSNR1_24998us_2025-04-15T111238_raw.hyspex";
hdr2 = "/Volumes/School/Thesis/Data/HySpex/cactusHalogen_002_VNIR_1800_SN00841_HSNR1_24998us_2025-04-15T112522_raw.hdr";
hyspex2 = "/Volumes/School/Thesis/Data/HySpex/cactusHalogen_002_VNIR_1800_SN00841_HSNR1_24998us_2025-04-15T112522_raw.hyspex";


%% Load HySpex radiance cube
function [cubeData, wavelengths] = loadHySpexCube(hdrFile, hyspexFile)
    info = enviinfo(hdrFile);
    cubeData = multibandread( ...
        hyspexFile, ...
        [info.Height, info.Width, info.Bands], ...
        info.DataType, ...
        info.HeaderOffset, ...
        info.Interleave, ...
        info.ByteOrder);                         
    wavelengths = info.Wavelength;
end

%% Load scans
[data1, wl] = loadHySpexCube(hdr1, hyspex1);
[data2, ~ ] = loadHySpexCube(hdr2, hyspex2);

%Select bands within 380–780 nm
bandMask = wl >= 380 & wl <= 780;
data1 = data1(:,:,bandMask);
data2 = data2(:,:,bandMask);
wl = wl(bandMask);  % Update wavelengths accordingly


%% Downsample both cubes 
scaleFactor = 0.5;
[origRows, origCols, B] = size(data1);

% Compute output sizes using round to match imresize’s ceil behavior 
dsRows = round(origRows * scaleFactor);
dsCols = round(origCols * scaleFactor);

% Preallocate downsampled cubes
data1_ds = zeros(dsRows, dsCols, B, 'like', data1);
data2_ds = zeros(dsRows, dsCols, B, 'like', data2);

for k = 1:B
    data1_ds(:,:,k) = imresize(data1(:,:,k), scaleFactor);    % outputs dsRows×dsCols 
    data2_ds(:,:,k) = imresize(data2(:,:,k), scaleFactor);
end


%% Estimate transformation using band 100
band_no = 110;
I1 = data1_ds(:,:,band_no);
I2 = data2_ds(:,:,band_no);

% Detect and extract features
points1 = detectSURFFeatures(I1);
points2 = detectSURFFeatures(I2);
[features1, validPoints1] = extractFeatures(I1, points1);
[features2, validPoints2] = extractFeatures(I2, points2);

% Match features
indexPairs = matchFeatures(features1, features2, 'Unique', true);
matchedPoints1 = validPoints1(indexPairs(:,1));
matchedPoints2 = validPoints2(indexPairs(:,2));

figure;
showMatchedFeatures(I1, I2, matchedPoints1, matchedPoints2, 'montage');
title('Matched SURF Features');

% Estimate transformation
[tform, ~, ~] = estimateGeometricTransform2D(...
    matchedPoints2, matchedPoints1, 'projective');

%%

registrationEstimator(I1, I2);                 
uiwait(warndlg(['Pick matching points, export transform as movingReg, then OK.'],'Setup'));

%%
tformObj = movingReg.Transformation;                          % geometric transform
Rfixed   = movingReg.SpatialRefObj;                           % spatial ref of fixed image

%%
tform = tformObj;
%% 3) Get the size of the images and their spatial alignment
[H1, W1] = size(I1);
[H2, W2] = size(I2);

% Calculate the overlap position based on the tform
% Use outputLimits to find the range of the images after transformation
[xlim1, ylim1] = outputLimits(projective2d(eye(3)), [1 W1], [1 H1]);
[xlim2, ylim2] = outputLimits(tform, [1 W2], [1 H2]);

% Find where the second image should start based on the overlap
xMin = floor(min([xlim1 xlim2])); 
xMax = ceil(max([xlim1 xlim2]));
yMin = floor(min([ylim1 ylim2]));
yMax = ceil(max([ylim1 ylim2]));

% Adjust the output size (width should be the sum of the two image widths minus overlap)
outputWidth = xMax - xMin;  % New width to avoid double overlap
outputHeight = yMax - yMin;  % Height remains the same for both images

% Define the output reference for the new canvas size
outputRef = imref2d([outputHeight, outputWidth], [xMin xMax], [yMin yMax]);

%% 4) Apply the transformations (warp the images)
I1w = imwarp(I1, affine2d(eye(3)), 'OutputView', outputRef);  % First image (no transformation)
I2w = imwarp(I2, tform, 'OutputView', outputRef);  % Second image with the estimated transform

%% 5) Overlay the images (avoiding overlap)
panorama = zeros(outputHeight, outputWidth, 'like', I1w);  % Initialize the canvas

% Place the first image in the left part of the canvas
panorama(1:size(I1w,1), 1:size(I1w,2)) = I1w;  

% Only place the second image where it has non-zero pixels
mask2 = I2w > 0;
panorama(mask2) = I2w(mask2);

%% 6) Display the result
figure;
imshow(panorama, []);
title(sprintf('Stitched Images', band_no));





