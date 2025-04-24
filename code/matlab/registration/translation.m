%% Full Translation-Based Hyperspectral Cube Stitching Script

clc; clear; close all;

%% 1) File paths
hdr1    = "/Volumes/School/Thesis/Data/HySpex/cactusHalogen_001_VNIR_1800_SN00841_HSNR1_24998us_2025-04-15T111238_raw.hdr";
hyspex1 = "/Volumes/School/Thesis/Data/HySpex/cactusHalogen_001_VNIR_1800_SN00841_HSNR1_24998us_2025-04-15T111238_raw.hyspex";
hdr2    = "/Volumes/School/Thesis/Data/HySpex/cactusHalogen_002_VNIR_1800_SN00841_HSNR1_24998us_2025-04-15T112522_raw.hdr";
hyspex2 = "/Volumes/School/Thesis/Data/HySpex/cactusHalogen_002_VNIR_1800_SN00841_HSNR1_24998us_2025-04-15T112522_raw.hyspex";

%% 2) Helper function: load a HySpex radiance cube
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

%% 3) Load and filter to visible bands (380–780 nm)
[data1, wl] = loadHySpexCube(hdr1, hyspex1);   % moving half
[data2, ~ ] = loadHySpexCube(hdr2, hyspex2);   % fixed half

bandMask = wl >= 380 & wl <= 780;
data1    = data1(:,:,bandMask);
data2    = data2(:,:,bandMask);
wl       = wl(bandMask);
[H, W, B] = size(data1);

%% 4) Downsample both halves by 50%
scaleFactor = 0.5;
H1 = round(H * scaleFactor);
W1 = round(W * scaleFactor);

data1_ds = zeros(H1, W1, B, 'like', data1);
data2_ds = zeros(H1, W1, B, 'like', data2);
for k = 1:B
    data1_ds(:,:,k) = imresize(data1(:,:,k), scaleFactor);
    data2_ds(:,:,k) = imresize(data2(:,:,k), scaleFactor);
end


%% 3) Find horizontal shift with normxcorr2
midBand = round(B / 2);
mov = im2single(data1_ds(:, :, midBand));
fix = im2single(data2_ds(:, :, midBand));

C = normxcorr2(mov, fix);
[~, idx] = max(C(:));
[ypeak, xpeak] = ind2sub(size(C), idx);
xoff = xpeak - size(mov, 2);  % Horizontal offset
fprintf('Horizontal offset: xoff = %d px\n', xoff);

%% 4) Crop and concatenate based on offset
absX = round(abs(xoff)/2) - 26;
if xoff < 0
    % Moving image is to the left of the fixed image
    fixedC = data2_ds(:, (absX):end, :);
    movingC = data1_ds;
    stitchedData = cat(2, movingC, fixedC);
elseif xoff > 0
    % Moving image is to the right of the fixed image
    fixedC = data2_ds;
    movingC = data1_ds(:, (1 + xoff):end, :);
    stitchedData = cat(2, fixedC, movingC);
else
    % No horizontal offset; images are aligned side by side
    stitchedData = cat(2, data2_ds, data1_ds);
end

%% 5) Create hypercube and visualize
stitchedCube = hypercube(stitchedData, wl);
rgbVis = colorize(stitchedCube, 'Method', 'RGB', 'ContrastStretching', true);
figure; imshow(rgbVis);
title('Correctly Ordered Stitched Hyperspectral Cube');




