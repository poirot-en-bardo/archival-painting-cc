%% File Paths
hdr1 = "/Volumes/School/Thesis/Data/HySpex/cactusHalogen_001_VNIR_1800_SN00841_HSNR1_24998us_2025-04-15T111238_raw.hdr";
hyspex1 = "/Volumes/School/Thesis/Data/HySpex/cactusHalogen_001_VNIR_1800_SN00841_HSNR1_24998us_2025-04-15T111238_raw.hyspex";
hdr2 = "/Volumes/School/Thesis/Data/HySpex/cactusHalogen_002_VNIR_1800_SN00841_HSNR1_24998us_2025-04-15T112522_raw.hdr";
hyspex2 = "/Volumes/School/Thesis/Data/HySpex/cactusHalogen_002_VNIR_1800_SN00841_HSNR1_24998us_2025-04-15T112522_raw.hyspex";

%% Load HySpex Radiance Cubes
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

[data1, wl] = loadHySpexCube(hdr1, hyspex1);
[data2, ~] = loadHySpexCube(hdr2, hyspex2);

%% Select Bands within 380–780 nm
bandMask = wl >= 380 & wl <= 780;
data1 = data1(:,:,bandMask);
data2 = data2(:,:,bandMask);
wl = wl(bandMask);  % Update wavelengths accordingly

%% Downsample Both Cubes
scaleFactor = 0.5;
[origRows, origCols, B] = size(data1);

dsRows = round(origRows * scaleFactor);
dsCols = round(origCols * scaleFactor);

data1_ds = zeros(dsRows, dsCols, B, 'like', data1);
data2_ds = zeros(dsRows, dsCols, B, 'like', data2);

for k = 1:B
    data1_ds(:,:,k) = imresize(data1(:,:,k), scaleFactor);    % outputs dsRows×dsCols 
    data2_ds(:,:,k) = imresize(data2(:,:,k), scaleFactor);
end

%% Choose a Band
band_no = 100;
I1 = data1_ds(:,:,band_no);
I2 = data2_ds(:,:,band_no);

%% Cross-Correlation to Find Horizontal Overlap
templateWidth = 50;
template = I1(:, end - templateWidth + 1:end); % last columns of I1

searchRegion = I2(:, 1:size(I1,2));  % first portion of I2

cc = normxcorr2(template, searchRegion);

[~, maxIdx] = max(cc(:));
[ypeak, xpeak] = ind2sub(size(cc), maxIdx);

matchX_in_searchRegion = xpeak - templateWidth + 1;

overlap = size(searchRegion,2) - matchX_in_searchRegion + 1;
% overlap = xpeak;
%% Compute Canvas Size
[H1, W1] = size(I1);
[H2, W2] = size(I2);

newWidth = W1 + W2 - overlap;
newHeight = max(H1, H2);

stitched = zeros(newHeight, newWidth, 'like', I1);

%% Place Full I1
stitched(1:H1, 1:W1) = I1;

%% Place Only Non-Overlapping Part of I2
stitched(1:H2, W1+1:end) = I2(:, overlap+1:end);

%% Show Result
figure;
imshow(stitched, []);
title('Stitched Image (Side-by-Side without Overlap)');