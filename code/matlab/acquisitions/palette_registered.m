clear; close all;
% cubeFile1 = '/home/oem/eliza/data/processed/reflectance/before/cactus_halogen_reflectance_full.hdr';
% cubeFile2 = '/home/oem/eliza/data/processed/reflectance/after/cactus_halogen_reflectance_after_full.hdr';
cubeFile1 = '/home/oem/eliza/data/processed/reflectance/before/yoda_halogen_reflectance_full.hdr';
cubeFile2 = '/home/oem/eliza/data/processed/reflectance/after/yoda_halogen_reflectance_after_full.hdr';

hcube1 = hypercube(cubeFile1); cube1 = hcube1.DataCube; wl1 = hcube1.Wavelength;
hcube2 = hypercube(cubeFile2); cube2 = hcube2.DataCube; wl2 = hcube2.Wavelength;
valid_idx1 = find(wl1 >= 380 & wl1 <= 780);
cube1 = cube1(:,:,valid_idx1); wl1 = wl1(valid_idx1);
valid_idx2 = find(wl2 >= 380 & wl2 <= 780);
cube2 = cube2(:,:,valid_idx2); wl2 = wl2(valid_idx2);

%%
% Crop 
figure; imagesc(cube1(:,:,10)); axis image; colormap gray; title('Crop FILM');
roi1 = drawrectangle; r1 = round(roi1.Position);
x1a = max(1, r1(1)); y1a = max(1, r1(2));
x1b = min(size(cube1,2), x1a + r1(3) - 1); y1b = min(size(cube1,1), y1a + r1(4) - 1);
cube1_crop = cube1(y1a:y1b, x1a:x1b, :);

figure; imagesc(cube2(:,:,10)); axis image; colormap gray; title('Crop PAINTING');
roi2 = drawrectangle; r2 = round(roi2.Position);
x2a = max(1, r2(1)); y2a = max(1, r2(2));
x2b = min(size(cube2,2), x2a + r2(3) - 1); y2b = min(size(cube2,1), y2a + r2(4) - 1);
cube2_crop = cube2(y2a:y2b, x2a:x2b, :);

% Subsample both to same size
scale = 0.2; % or whatever you want
sz = [round(size(cube1_crop,1)*scale), round(size(cube1_crop,2)*scale), size(cube1_crop,3)];
cube1_sub = imresize3(cube1_crop, sz);
cube2_sub = imresize3(cube2_crop, sz);

%% registration

bandToShow = 85; 

I1 = mat2gray(cube1_sub(:,:,bandToShow)); % Reference (before)
I2 = mat2gray(cube2_sub(:,:,bandToShow)); % To register (after)

% --- Detect features
points1 = detectSURFFeatures(I1);
points2 = detectSURFFeatures(I2);

% --- Extract feature descriptors
[features1, valid_points1] = extractFeatures(I1, points1);
[features2, valid_points2] = extractFeatures(I2, points2);

% --- Match features
indexPairs = matchFeatures(features1, features2, 'MaxRatio', 0.5, 'Unique', true);
matchedPoints1 = valid_points1(indexPairs(:,1));
matchedPoints2 = valid_points2(indexPairs(:,2));

% --- Estimate transformation (affine or similarity)
[tform, inlierIdx] = estimateGeometricTransform2D(matchedPoints2, matchedPoints1, 'affine', ...
    'MaxNumTrials', 5000, 'MaxDistance', 3);
% You can try 'similarity' or 'projective' too

% --- Visualize matches (optional)
figure; showMatchedFeatures(I1, I2, matchedPoints1, matchedPoints2, 'montage');
title('Matched SURF points');

% --- Register "after" cube (all bands)
Rfixed = imref2d(size(I1));
cube2_reg = zeros(size(cube1_sub));
for b = 1:size(cube2_sub,3)
    cube2_reg(:,:,b) = imwarp(cube2_sub(:,:,b), tform, 'OutputView', Rfixed, 'Interp', 'bilinear');
end

% --- Overlap mask (valid area in both cubes)
valid1 = any(cube1_sub > 0.01, 3);
valid2 = any(cube2_reg > 0.01, 3);
valid = valid1 & valid2;

[rows, cols] = find(valid);
rowmin = min(rows); rowmax = max(rows);
colmin = min(cols); colmax = max(cols);

cube1_final = cube1_sub(rowmin:rowmax, colmin:colmax, :);
cube2_final = cube2_reg(rowmin:rowmax, colmin:colmax, :);

% Optional: show overlays to check alignment
figure; imshowpair(mat2gray(cube1_final(:,:,bandToShow)), mat2gray(cube2_final(:,:,bandToShow)));
title('After Automatic Registration');

%% palette extraction
% Divide into blocks
[H, W, B] = size(cube1_final);
block_h = 10; block_w = 10;
n_blocks_y = floor(H/block_h);
n_blocks_x = floor(W/block_w);

meanSpectra1 = zeros(n_blocks_y*n_blocks_x, B);
meanSpectra2 = zeros(n_blocks_y*n_blocks_x, B);

idx = 1;
for i = 1:n_blocks_y
    for j = 1:n_blocks_x
        r_idx = (i-1)*block_h+1:i*block_h;
        c_idx = (j-1)*block_w+1:j*block_w;

        block1 = reshape(cube1_final(r_idx, c_idx, :), [], B);
        block2 = reshape(cube2_final(r_idx, c_idx, :), [], B);

        meanSpectra1(idx,:) = mean(block1, 1, 'omitnan');
        meanSpectra2(idx,:) = mean(block2, 1, 'omitnan');
        idx = idx + 1;
    end
end


%%
% Paths to CMF and D65 files
cmf = importdata('../../../data/CIE2degCMFs_1931.txt');
ill = importdata('../../../data/CIE_D65.txt');
cmf_wl = cmf(:,1);
cmf_xyz = cmf(:,2:4);
ill_wl = ill(:,1);
ill_val = ill(:,2);
% Interpolate illuminant and CMFs to your wavelengths
ill_interp = interp1(ill_wl, ill_val, wl1, 'linear', 'extrap');
cmf_interp = interp1(cmf_wl, cmf_xyz, wl1, 'linear', 'extrap');

nColors = size(meanSpectra1, 1);

% Find white point for scaling XYZ to Lab (Y=100 for white)
white_refl = ones(1, numel(wl1));
white_xyz = (white_refl .* ill_interp(:)') * cmf_interp;
white_xyz = white_xyz / sum(ill_interp(:)' .* cmf_interp(:,2)');
ill_interp = ill_interp(:)';
% Convert spectra to XYZ, then scale, then to Lab & RGB
XYZs1 = (meanSpectra1 .* ill_interp) * cmf_interp;
XYZs2 = (meanSpectra2 .* ill_interp) * cmf_interp;
XYZs1 = XYZs1 ./ sum(ill_interp .* cmf_interp(:,2)', 2); % Normalize
XYZs2 = XYZs2 ./ sum(ill_interp .* cmf_interp(:,2)', 2);

XYZs1_scaled = XYZs1 * (100 / white_xyz(2));
XYZs2_scaled = XYZs2 * (100 / white_xyz(2));

Lab1 = xyz2lab(XYZs1_scaled, 'WhitePoint', 'd65');
Lab2 = xyz2lab(XYZs2_scaled, 'WhitePoint', 'd65');

RGBs1 = max(0, min(1, xyz2rgb(XYZs1, 'ColorSpace','srgb')));
RGBs2 = max(0, min(1, xyz2rgb(XYZs2, 'ColorSpace','srgb')));

% Block palette grid size
grid_h = n_blocks_y;
grid_w = n_blocks_x;
patchSize = 40;

img1 = ones(grid_h * patchSize, grid_w * patchSize, 3);
img2 = ones(grid_h * patchSize, grid_w * patchSize, 3);

for k = 1:nColors
    row = floor((k-1)/grid_w);
    col = mod((k-1), grid_w);
    r_idx = (row*patchSize+1):((row+1)*patchSize);
    c_idx = (col*patchSize+1):((col+1)*patchSize);
    img1(r_idx, c_idx, :) = repmat(reshape(RGBs1(k,:),1,1,3), patchSize, patchSize, 1);
    img2(r_idx, c_idx, :) = repmat(reshape(RGBs2(k,:),1,1,3), patchSize, patchSize, 1);
end


%
figure; imshow(img1); title('Palette from Cube 1 (Before)');
figure; imshow(img2); title('Palette from Cube 2 (After)');

%% k means
% After registration and masking to valid overlap:
% cube1_final, cube2_final: [H x W x B]
pixels1 = reshape(cube1_final, [], size(cube1_final,3)); % [N x B]
pixels2 = reshape(cube2_final, [], size(cube2_final,3)); % [N x B]
valid = all(~isnan(pixels1),2) & all(~isnan(pixels2),2); % Only fully valid pixels
pixels1 = pixels1(valid, :);
pixels2 = pixels2(valid, :);

nColors = 60;
[cluster_idx, ~] = kmeans(pixels1, nColors, 'Replicates', 5, 'MaxIter', 1000);

meanSpectra1 = zeros(nColors, size(pixels1,2));
meanSpectra2 = zeros(nColors, size(pixels2,2));
for k = 1:nColors
    meanSpectra1(k,:) = mean(pixels1(cluster_idx==k,:), 1, 'omitnan');
    meanSpectra2(k,:) = mean(pixels2(cluster_idx==k,:), 1, 'omitnan');
end

%% Step 7: Save
save('palette_yoda_before1.mat', 'meanSpectra1', 'wl1');
save('palette_yoda_after1.mat', 'meanSpectra2', 'wl2');

%%

cmf = importdata('../../../data/CIE2degCMFs_1931.txt');  % or use your path
ill = importdata('../../../data/CIE_D65.txt');
cmf_wl = cmf(:,1);
cmf_xyz = cmf(:,2:4);
ill_wl = ill(:,1);
ill_val = ill(:,2);

nColors = size(meanSpectra1, 1);

% Interpolate illuminant and cmfs to match your wavelengths
ill_interp = interp1(ill_wl, ill_val, wl1, 'linear', 'extrap');
cmf_interp = interp1(cmf_wl, cmf_xyz, wl1, 'linear', 'extrap');

% Function to convert spectra to xyz then to rgb
spectra2rgb = @(refl) xyz2rgb((refl .* ill_interp') * cmf_interp ./ sum(ill_interp' .* cmf_interp(:,2)), 'ColorSpace','srgb');

% Convert each mean spectrum to RGB
RGBs1 = zeros(nColors, 3);
RGBs2 = zeros(nColors, 3);
for k = 1:nColors
    refl = meanSpectra1(k,:); % 1 x N
    weighted = refl .* ill_interp(:)'; % 1 x N
    xyz1 = weighted * cmf_interp; % 1 x 3
    xyz1 = xyz1 / sum(ill_interp(:)' .* cmf_interp(:,2)');
    xyz1 = max(0, xyz1);
    RGBs1(k,:) = xyz2rgb(xyz1, 'ColorSpace','srgb');

    refl = meanSpectra2(k,:);
    weighted = refl .* ill_interp(:)';
    xyz2 = weighted * cmf_interp;
    xyz2 = xyz2 / sum(ill_interp(:)' .* cmf_interp(:,2)');
    xyz2 = max(0, xyz2);
    RGBs2(k,:) = xyz2rgb(xyz2, 'ColorSpace','srgb');
end


% Clamp RGBs to [0 1]
RGBs1 = max(0, min(1, RGBs1));
RGBs2 = max(0, min(1, RGBs2));

% --- (2) Make a grid image for each palette


nColors = 60;
grid_w = 10;
grid_h = 6;
patchSize = 40;

% --- img1 and img2 initialised as white ---
img1 = ones(grid_h * patchSize, grid_w * patchSize, 3); % Cube 1
img2 = ones(grid_h * patchSize, grid_w * patchSize, 3); % Cube 2 (matched!)

for k = 1:nColors
    row = floor((k-1)/grid_w);
    col = mod((k-1), grid_w);
    r_idx = (row*patchSize+1):((row+1)*patchSize);
    c_idx = (col*patchSize+1):((col+1)*patchSize);
    patch_rgb1 = repmat(reshape(RGBs1(k,:), 1, 1, 3), patchSize, patchSize, 1);
    patch_rgb2 = repmat(reshape(RGBs2(k,:), 1, 1, 3), patchSize, patchSize, 1);
    img1(r_idx, c_idx, :) = patch_rgb1;               % before
    img2(r_idx, c_idx, :) = patch_rgb2;               % after (matched!)
end




%%
figure; imshow(img1); title('Palette Yoda (Before)');
figure; imshow(img2); title('Palette Yoda (After)');


%%
% After you have Lab1 and Lab2 as [numBlocks x 3] or [nColors x 3] arrays
dE = deltaE2000(Lab1, Lab2);

% Put dE into a grid
dE_grid = reshape(dE, grid_h, grid_w);
figure; imagesc(dE_grid);
axis image off;
colormap(jet(255)); colorbar; clim([0 10]);
title('\DeltaE_{00} between palette patches');
for k = 1:numel(dE)
    [row, col] = ind2sub([grid_h, grid_w], k);
    text(col, row, sprintf('%.1f', dE(k)), ...
        'HorizontalAlignment','center','Color','w','FontSize',10,'FontWeight','bold');
end


%%
%saving registered cubes
outFolder = '/home/oem/eliza/data/processed/reflectance/registered';
basename1 = 'yoda_reg_before1';
basename2 = 'yoda_reg_after1';

datFile1 = fullfile(outFolder, [basename1 '.dat']);
hdrFile1 = fullfile(outFolder, [basename1 '.hdr']);
datFile2 = fullfile(outFolder, [basename2 '.dat']);
hdrFile2 = fullfile(outFolder, [basename2 '.hdr']);

% Save the cubes as .dat (BIL)
multibandwrite(single(cube1_final), datFile1, 'bil');
multibandwrite(single(cube2_final), datFile2, 'bil');

% Prepare header info
[rows, cols, bands] = size(cube1_final);
hdr.samples = cols;
hdr.lines = rows;
hdr.bands = bands;
hdr.header_offset = 0;
hdr.data_type = 4;        % 4 = single (float)
hdr.interleave = 'bil';
hdr.byte_order = 0;       % 0 = little-endian
hdr.wavelength = wl1(:)'; % Or wl2(:)' for cube2

write_envi_hdr(hdrFile1, hdr);
write_envi_hdr(hdrFile2, hdr);


function write_envi_hdr(hdrfile, hdr)
    fid = fopen(hdrfile, 'w');
    fprintf(fid, 'ENVI\n');
    fprintf(fid, 'samples = %d\n', hdr.samples);
    fprintf(fid, 'lines   = %d\n', hdr.lines);
    fprintf(fid, 'bands   = %d\n', hdr.bands);
    fprintf(fid, 'header offset = %d\n', hdr.header_offset);
    fprintf(fid, 'file type = {ENVI Standard}\n');
    fprintf(fid, 'data type = %d\n', hdr.data_type);
    fprintf(fid, 'byte order = %d\n', hdr.byte_order); 
    fprintf(fid, 'interleave = %s\n', hdr.interleave);
    if isfield(hdr, 'wavelength')
        fprintf(fid, 'wavelength = { ');
        fprintf(fid, '%g, ', hdr.wavelength(1:end-1));
        fprintf(fid, '%g }\n', hdr.wavelength(end));
    end
    fclose(fid);
end
