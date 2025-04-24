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


%% Interactive registration on mid-spectrum band
midBand   = round(B/2);
midBand   = 100;
movingImg = data1_ds(:,:,midBand);
fixedImg  = data2_ds(:,:,midBand);

registrationEstimator(movingImg, fixedImg);                 
uiwait(warndlg(['Pick matching points, export transform as movingReg, then OK.'],'Setup'));

tformObj = movingReg.Transformation;                          % geometric transform
Rfixed   = movingReg.SpatialRefObj;                           % spatial ref of fixed image

%% Compute unified canvas limits
xLimits = [1, size(data1_ds,2)];
yLimits = [1, size(data1_ds,1)];
[xW, yW] = outputLimits(tformObj, xLimits, yLimits);
xF       = Rfixed.XWorldLimits;
yF       = Rfixed.YWorldLimits;

xMin = min(xF(1), xW(1));  xMax = max(xF(2), xW(2));
yMin = min(yF(1), yW(1));  yMax = max(yF(2), yW(2));
outW = round(xMax - xMin);
outH = round(yMax - yMin);

outputView = imref2d([outH, outW], [xMin xMax], [yMin yMax]);

%% 7) Warp ONLY the moving downsampled cube; pad the fixed one
warpedMoving = imwarp(data1_ds, tformObj, 'OutputView', outputView);
warpedFixed  = imwarp(data2_ds, affine2d(eye(3)), 'OutputView', outputView);

%% 8) Merge: fill gaps in moving with fixed, then blend overlap
maskNoData = (warpedMoving == 0);
stitched   = warpedMoving;
stitched(maskNoData) = warpedFixed(maskNoData);

ov = ~maskNoData & (warpedFixed ~= 0);
stitched(ov) = 0.5 * (double(warpedMoving(ov)) + double(warpedFixed(ov)));

%% 9) Build the final full-painting hypercube from downsampled data
stitchedCube = hypercube(stitched, wl);

% 10) (Optional) Visualize an RGB composite of the stitched result
rgbVis = colorize(stitchedCube, 'Method','RGB','ContrastStretching',true);
figure; imshow(rgbVis);
title('Full-Painting Hyperspectral Cube (Downsampled & Stitched)');

