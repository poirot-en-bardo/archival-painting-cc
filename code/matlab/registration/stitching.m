
hdr1 = "/Volumes/School/Thesis/Data/HySpex/cactusHalogen_001_VNIR_1800_SN00841_HSNR1_24998us_2025-04-15T111238_raw.hdr";
hyspex1 = "/Volumes/School/Thesis/Data/HySpex/cactusHalogen_001_VNIR_1800_SN00841_HSNR1_24998us_2025-04-15T111238_raw.hyspex";
hdr2 = "/Volumes/School/Thesis/Data/HySpex/cactusHalogen_002_VNIR_1800_SN00841_HSNR1_24998us_2025-04-15T112522_raw.hdr";
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

%% 3) Load both scans
[data1, wl] = loadHySpexCube(hdr1, hyspex1);
[data2, ~ ] = loadHySpexCube(hdr2, hyspex2);

%% 4) Interactive registration on a visible band
midBand   = round(numel(wl)/2);
movingImg = data1(:,:,midBand);
fixedImg  = data2(:,:,midBand);

registrationEstimator(movingImg, fixedImg);
uiwait(warndlg( ...
    'Select matching points in the overlap, export transform as ''movingReg'', then click OK.', ...
    'Stitch Setup'));

tformObj = movingReg.Transformation;    % projective2d object
Rfixed   = movingReg.SpatialRefObj;     % imref2d for scan2

%% 5) Compute combined output limits
[xInMin, xInMax] = deal([1, size(data1,2)]);
[yInMin, yInMax] = deal([1, size(data1,1)]);
[xWarped, yWarped] = outputLimits(tformObj, xInMin, yInMax);
xFixed = Rfixed.XWorldLimits;
yFixed = Rfixed.YWorldLimits;

xMin = min(xFixed(1), xWarped(1));
xMax = max(xFixed(2), xWarped(2));
yMin = min(yFixed(1), yWarped(1));
yMax = max(yFixed(2), yWarped(2));

width  = round(xMax - xMin);
height = round(yMax - yMin);

outputView = imref2d([height, width], [xMin xMax], [yMin yMax]);

%% 6) Warp and stitch using block-wise processing
B = size(data1,3);  % Number of spectral bands
blockSize = 10;     % Number of bands to process at a time

identityTform = affine2d(eye(3));

% Preallocate a cell array to store stitched blocks
stitchedBlocks = cell(1, ceil(B / blockSize));

for i = 1:blockSize:B
    idx = i:min(i+blockSize-1, B);
    
    % Warp blocks
    warped1 = zeros(height, width, numel(idx), 'like', data1);
    warped2 = zeros(height, width, numel(idx), 'like', data2);
    
    for k = 1:numel(idx)
        warped1(:,:,k) = imwarp(data1(:,:,idx(k)), tformObj, ...
                                'OutputView', outputView);
        warped2(:,:,k) = imwarp(data2(:,:,idx(k)), identityTform, ...
                                'OutputView', outputView);
    end
    
    % Stitch blocks
    maskNoData = (warped1 == 0);
    blockStitched = warped1;
    blockStitched(maskNoData) = warped2(maskNoData);
    
    overlap = ~maskNoData & (warped2 ~= 0);
    blockStitched(overlap) = 0.5 * (double(warped1(overlap)) + double(warped2(overlap)));
    
    % Store the stitched block
    stitchedBlocks{ceil(i / blockSize)} = blockStitched;
end

% Concatenate all stitched blocks along the third dimension
stitched = cat(3, stitchedBlocks{:});

%% 7) Build the output hypercube
stitchedCube = hypercube(stitched, wl);

%% 8) (Optional) Visualize an RGB composite
rgbVis = colorize(stitchedCube, ...
    'Method', 'RGB', ...
    'ContrastStretching', true);
figure; imshow(rgbVis);
title('Stitched Hyperspectral Painting (RGB Composite)');