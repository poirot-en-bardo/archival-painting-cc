

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


%% Use band 100 to estimate transform
band_no   = 100;
gray1 = data1_ds(:,:,band_no);
gray2 = data2_ds(:,:,band_no);


% Ensure images are in double precision for processing
I1 = data1_ds(:,:,band_no);
I2 = data2_ds(:,:,band_no);

% Estimate overlap using normalized cross-correlation
c = normxcorr2(I2, I1);
[~, imax] = max(abs(c(:)));
[ypeak, xpeak] = ind2sub(size(c), imax);

% Calculate the offset between the images
yoffset = ypeak - size(I2, 1);
xoffset = xpeak - size(I2, 2);

% Determine the size of the stitched image
rows = max(size(I1, 1), yoffset + size(I2, 1));
cols = max(size(I1, 2), xoffset + size(I2, 2));
stitched_image = zeros(rows, cols);

% Place I1 in the stitched image
stitched_image(1:size(I1, 1), 1:size(I1, 2)) = I1;



% Overlay I2 onto the stitched image with blending
for y = 1:size(I2, 1)
    for x = 1:size(I2, 2)
        row = yoffset + y;
        col = xoffset + x;
       
        % Directly assign non-overlapping pixels
        stitched_image(row, col) = I2(y, x);
        mask(row, col) = 1;
        
    end
end

% Display the stitched image
imshow(stitched_image, []);
title('Stitched Image');
