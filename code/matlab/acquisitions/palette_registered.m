clear; close all;
% cubeFile1 = '/home/oem/eliza/data/processed/reflectance/before/cactus_halogen_reflectance_full.hdr';
% cubeFile2 = '/home/oem/eliza/data/processed/reflectance/after/cactus_halogen_reflectance_after_full.hdr';
% cubeFile1 = '/home/oem/eliza/data/processed/reflectance/before/yoda_halogen_reflectance_full.hdr';
% cubeFile2 = '/home/oem/eliza/data/processed/reflectance/after/yoda_halogen_reflectance_after_full.hdr';
cubeFile1 = '/home/oem/eliza/data/reflectance/registered/cactus_reflectance_before.mat';
cubeFile2 = '/home/oem/eliza/data/reflectance/registered/cactus_reflectance_after_reg.mat';
%%
hcube_before = load(cubeFile1);
hcube_after = load(cubeFile2);
%%
hcube1 = hypercube(cubeFile1); cube1 = hcube1.DataCube; wl1 = hcube1.Wavelength;
hcube2 = hypercube(cubeFile2); cube2 = hcube2.DataCube; wl2 = hcube2.Wavelength;
valid_idx1 = find(wl1 >= 380 & wl1 <= 780);
cube1 = cube1(:,:,valid_idx1); wl1 = wl1(valid_idx1);
valid_idx2 = find(wl2 >= 380 & wl2 <= 780);
cube2 = cube2(:,:,valid_idx2); wl2 = wl2(valid_idx2);

%%
% Crop 
% For cube 1
figure; imagesc(cube1(:,:,80)); axis image; colormap gray; title('Crop FILM');
roi1 = drawrectangle('InteractionsAllowed','all');
% Wait for you to double-click inside rectangle to confirm
wait(roi1);
r1 = round(roi1.Position);
x1a = max(1, r1(1)); y1a = max(1, r1(2));
x1b = min(size(cube1,2), x1a + r1(3) - 1); y1b = min(size(cube1,1), y1a + r1(4) - 1);
cube1_crop = cube1(y1a:y1b, x1a:x1b, :);

% For cube 2
figure; imagesc(cube2(:,:,80)); axis image; colormap gray; title('Crop PAINTING');
roi2 = drawrectangle('InteractionsAllowed','all');
wait(roi2);
r2 = round(roi2.Position);
x2a = max(1, r2(1)); y2a = max(1, r2(2));
x2b = min(size(cube2,2), x2a + r2(3) - 1); y2b = min(size(cube2,1), y2a + r2(4) - 1);
cube2_crop = cube2(y2a:y2b, x2a:x2b, :);

%%
% Subsample both to same size
scale = 0.5; % or whatever you want
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
indexPairs = matchFeatures(features1, features2, 'MaxRatio', 0.3, 'Unique', true);
matchedPoints1 = valid_points1(indexPairs(:,1));
matchedPoints2 = valid_points2(indexPairs(:,2));

% --- Estimate transformation (affine or similarity)
[tform, inlierIdx] = estimateGeometricTransform2D(matchedPoints2, matchedPoints1, 'affine', ...
    'MaxNumTrials', 5000, 'MaxDistance', 3);

% --- Visualize matches 
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
%
% Optional: show overlays to check alignment
figure; imshowpair(mat2gray(cube1_final(:,:,bandToShow)), mat2gray(cube2_final(:,:,bandToShow)));
title('After Automatic Registration');



%% using already-registered cubes 

% img_path1 = '/home/oem/eliza/data/processed/reflectance/registered/small/yoda_reg_before1.hdr';
img_path1 = '/home/oem/eliza/data/processed/reflectance/registered/cactus_reg_before1.hdr';
hcube1= hypercube(img_path1);
hcube2 = hypercube('/home/oem/eliza/data/processed/reflectance/registered/cactus_reg_after1.hdr');
% hcube1 = hypercube(img_path1);
% hcube2 = hypercube('/home/oem/eliza/data/processed/reflectance/registered/small/yoda_reg_after1.hdr');

cube1_final = hcube1.DataCube;
cube2_final = hcube2.DataCube;
wl1 = hcube1.Wavelength;


%%
% cube1_final = fliplr(cube1_final);
% cube2_final = fliplr(cube2_final);


%% K means only on non-specular pixels

mask1 = ~all(cube1_final > 0.98, 3);
mask2 = ~all(cube2_final > 0.98, 3);

%  Make intersection mask of valid pixels
min_h = min(size(cube1_final,1), size(cube1_final,1));
min_w = min(size(cube2_final,2), size(cube2_final,2));
cube1_final = cube1_final(1:min_h,1:min_w,:);
cube2_final = cube2_final(1:min_h,1:min_w,:);
mask1 = mask1(1:min_h,1:min_w);
mask2 = mask2(1:min_h,1:min_w);

mask = mask1 & mask2; % Intersection

% extracting pixels
pixels1 = reshape(cube1_final, [], size(cube1_final,3));
pixels2 = reshape(cube2_final, [], size(cube2_final,3));
mask_flat = mask(:);
pixels1_masked = pixels1(mask_flat, :);
pixels2_masked = pixels2(mask_flat, :);

rng(42);
nColors = 80;
[cluster_idx, ~] = kmeans(pixels1_masked, nColors, 'Replicates', 5, 'MaxIter', 1000);

meanSpectra1 = zeros(nColors, size(pixels1_masked,2));
meanSpectra2 = zeros(nColors, size(pixels2,2));
for k = 1:nColors
    meanSpectra1(k,:) = mean(pixels1_masked(cluster_idx==k,:), 1, 'omitnan');
    meanSpectra2(k,:) = mean(pixels2_masked(cluster_idx==k,:), 1, 'omitnan');
end




%%

save_folder = '/home/oem/eliza/mac-shared/palette'; 
[~, img_name, ~] = fileparts(img_path1);

% Load CIE and D65 data
cmf = importdata('../../../data/CIE2degCMFs_1931.txt');
ill = importdata('../../../data/CIE_D50.txt');
cmf_interp = interp1(cmf(:,1), cmf(:,2:4), wl1, 'linear', 'extrap');
ill_interp = interp1(ill(:,1), ill(:,2), wl1, 'linear', 'extrap')';

nColors = size(meanSpectra1, 1);


% Calculate XYZ for all pixels
XYZ1 = ref2xyz(ill_interp(:), cmf_interp, meanSpectra1);  % returns [numPixels x 3]
XYZ2 = ref2xyz(ill_interp(:), cmf_interp, meanSpectra2);

% XYZ to Lab
Lab1 = xyz2lab(XYZ1);
Lab2 = xyz2lab(XYZ2);

XYZ1_norm = XYZ1 ./ 100;
XYZ2_norm = XYZ2 ./ 100;

RGBs1 = max(0, min(1, xyz2prophoto(XYZ1_norm, true)));
RGBs2 = max(0, min(1, xyz2prophoto(XYZ2_norm, true)));

Lab1_unsorted = Lab1;
Lab2_unsorted = Lab2;
RGBs1_unsorted = RGBs1;
RGBs2_unsorted = RGBs2;

%

% Compute chroma from Lab2 (use after, but you can use before if you prefer)
Chroma = sqrt(Lab1(:,2).^2 + Lab1(:,3).^2);
achromaticThresh = 10;   % adjust
isAchromatic = Chroma < achromaticThresh;

achromaticIdx = find(isAchromatic);
chromaticIdx = find(~isAchromatic);

% Cluster only chromatic (colored) patches
nGroups = 3;
if ~isempty(chromaticIdx)
    [grp, C] = kmeans(Lab1(chromaticIdx,2:3), nGroups, 'Replicates', 10);
else
    grp = [];
end

% Order: (1) all achromatics, sorted by lightness, (2) then colored by group/lightness
finalOrder = [];

% 1. Achromatics, by lightness
[~, achroOrd] = sort(Lab1(achromaticIdx,1), 'descend');
finalOrder = [finalOrder; achromaticIdx(achroOrd)];

% 2. Each chromatic group, by lightness
for g = 1:nGroups
    idx = chromaticIdx(grp==g);   % indices in original array
    [~, ord] = sort(Lab2(idx,1), 'descend');
    finalOrder = [finalOrder; idx(ord)];
end

%

%
% Now sort for display
Lab1 = Lab1_unsorted(finalOrder,:);
Lab2 = Lab2_unsorted(finalOrder,:);
RGBs1 = RGBs1_unsorted(finalOrder,:);
RGBs2 = RGBs2_unsorted(finalOrder,:);

%
% Palette grid display
patchSize = 40;
grid_w = 10; % columns
grid_h = nColors/grid_w;
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
figure; imshow(img1); title('Palette Yoda (Before)');
%
figure; imshow(img2); title('Palette Yoda (After)');


%%
% Delta E grid
% Delta E grid (keep this as row-wise, as palette)
dE = deltaE2000(Lab1, Lab2);
dE_grid = zeros(grid_h, grid_w);
for k = 1:nColors
    row = floor((k-1)/grid_w) + 1;
    col = mod((k-1), grid_w) + 1;
    dE_grid(row, col) = dE(k);
end

figure;
imagesc(dE_grid);
axis image off;
colormap(jet(255));
colorbar;
clim([0 10]);
title('\DeltaE_{00} between palette patches');
for k = 1:nColors
    row = floor((k-1)/grid_w) + 1;
    col = mod((k-1), grid_w) + 1;
    text(col, row, sprintf('%.1f', dE(k)), ...
        'HorizontalAlignment','center', 'Color','w', 'FontSize',10, 'FontWeight','bold');
end

save_path = fullfile(save_folder, ['deltaE_grid_' img_name '.png']);
saveas(gcf, save_path);

%% Save
% meanSpectra1 = meanSpectra1(finalOrder, :); %sorting
% meanSpectra2 = meanSpectra2(finalOrder, :);
% save('palette/palette_cactus_before1.mat', 'meanSpectra1', 'wl1');
% save('palette/palette_cactus_after1.mat', 'meanSpectra2', 'wl2');

%%
mean_spectra_before = meanSpectra1(finalOrder, :);
mean_spectra_after  = meanSpectra2(finalOrder, :);

xyz_before = XYZ1(finalOrder, :);
xyz_after  = XYZ2(finalOrder, :);

Lab_before = xyz2lab(xyz_before);
Lab_after  = xyz2lab(xyz_after);

% ΔE between before and after for each patch
deltaE = deltaE2000(Lab_before, Lab_after);

%%
save_path = fullfile(save_folder, ['palette_' img_name '.mat']);


save(save_path, ...
    'mean_spectra_before', ...
    'mean_spectra_after', ...
    'xyz_before', ...
    'xyz_after', ...
    'deltaE', ...
    '-v7.3');

fprintf('Saved palette patch data to: %s\n', save_path);

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

%%
% Display indices on palette (Before)
figure; imshow(img1); title('Palette (Before)');
hold on;
for k = 1:nColors
    row = floor((k-1)/grid_w);
    col = mod((k-1), grid_w);
    xpos = col*patchSize + patchSize/2;
    ypos = row*patchSize + patchSize/2;
    text(xpos, ypos, num2str(k), 'Color', 'w', 'FontSize', 14, ...
        'HorizontalAlignment','center', 'FontWeight','bold');
end
hold off;
%%
% Display indices on palette (After)
figure; imshow(img2); title('Palette (After)');
hold on;
for k = 1:nColors
    row = floor((k-1)/grid_w);
    col = mod((k-1), grid_w);
    xpos = col*patchSize + patchSize/2;
    ypos = row*patchSize + patchSize/2;
    text(xpos, ypos, num2str(k), 'Color', 'w', 'FontSize', 14, ...
        'HorizontalAlignment','center', 'FontWeight','bold');
end
hold off;

%%
% Overlay numbers in the same order as palette
figure; imagesc(dE_grid);
axis image off; colormap(jet(255)); colorbar; clim([0 10]);
title('\DeltaE_{00} between palette patches');
for k = 1:nColors
    row = floor((k-1)/grid_w) + 1;
    col = mod((k-1), grid_w) + 1;
    text(col, row, num2str(k), ...
        'HorizontalAlignment','center', 'Color','w', 'FontSize',14, 'FontWeight','bold');
end


%%
% --- Palette grid display with diagonal split: before (upper left), after (lower right) ---

patchSize = 40;
grid_w = 10;
grid_h = nColors / grid_w;
img_diag = ones(grid_h * patchSize, grid_w * patchSize, 3);

for k = 1:nColors
    row = floor((k-1)/grid_w);
    col = mod((k-1), grid_w);
    r_idx = (row*patchSize+1):((row+1)*patchSize);
    c_idx = (col*patchSize+1):((col+1)*patchSize);
    
    patch_before = repmat(reshape(RGBs1(k,:),1,1,3), patchSize, patchSize, 1);
    patch_after  = repmat(reshape(RGBs2(k,:),1,1,3), patchSize, patchSize, 1);
    
    % Create diagonal mask: upper left (before), lower right (after)
    [X, Y] = meshgrid(1:patchSize, 1:patchSize);
    diag_mask = Y < X;      % Lower triangle: after, Upper triangle (including diagonal): before

    patch = patch_before;   % Start with 'before'
    for c = 1:3
        temp = patch(:,:,c);
        temp2 = patch_after(:,:,c);   % Store patch_after channel in a temp variable
        temp(diag_mask) = temp2(diag_mask);
        patch(:,:,c) = temp;
    end
    
    img_diag(r_idx, c_idx, :) = patch;
end

figure;
imshow(img_diag);
title('Palette: Before/After');

%%
% Scale to uint16 (16-bit) - assuming your images are [0,1] double
img1_uint16 = uint16(round(img1 * 65535));
img2_uint16 = uint16(round(img2 * 65535));
img_diag_uint16 = uint16(round(img_diag * 65535));

% Specify your save folder and file names
saveProPhotoTIFF(img1_uint16,   fullfile(save_folder, [img_name '_palette_before.tif']));
saveProPhotoTIFF(img2_uint16,   fullfile(save_folder, [img_name '_palette_after.tif']));
saveProPhotoTIFF(img_diag_uint16, fullfile(save_folder, [img_name '_palette_diag.tif']));


%%
figure; imshow(img1); title('Palette Yoda (Before)');
hold on;
for k = 1:nColors
    row = floor((k-1)/grid_w);
    col = mod((k-1), grid_w);
    xpos = col*patchSize + patchSize/2;
    ypos = row*patchSize + patchSize/2;
    text(xpos, ypos, num2str(k), 'Color', 'w', 'FontSize', 14, ...
        'HorizontalAlignment','center', 'FontWeight','bold');
end
hold off;

figure;
imagesc(dE_grid);
axis image off; colormap(jet(255)); colorbar; clim([0 10]);
title('\DeltaE_{00} between palette patches');
for k = 1:nColors
    row = floor((k-1)/grid_w) + 1;
    col = mod((k-1), grid_w) + 1;
    text(col, row, num2str(k), ...
        'HorizontalAlignment','center', 'Color','w', 'FontSize',14, 'FontWeight','bold');
end


%% block-wise approach
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



% Calculate XYZ for all pixels
XYZ1 = ref2xyz(ill_interp(:), cmf_interp, meanSpectra1);  % returns [numPixels x 3]
XYZ2 = ref2xyz(ill_interp(:), cmf_interp, meanSpectra2);

% XYZ to Lab
Lab1 = xyz2lab(XYZ1);
Lab2 = xyz2lab(XYZ2);

% RGBs1 = max(0, min(1, xyz2rgb(XYZs1, 'ColorSpace','srgb')));
% RGBs2 = max(0, min(1, xyz2rgb(XYZs2, 'ColorSpace','srgb')));

RGBs1 = max(0, min(1, xyz2rgb(XYZ1, 'ColorSpace','srgb')));
RGBs2 = max(0, min(1, xyz2rgb(XYZ2, 'ColorSpace','srgb')));

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

%%
%Checking

data = load('/home/oem/eliza/mac-shared/palette/palette_cactus_reg_before1.mat');
RGB_before = max(0, min(1, xyz2rgb(data.xyz_before ./ 100, 'ColorSpace', 'srgb')));
RGB_after  = max(0, min(1, xyz2rgb(data.xyz_after  ./ 100, 'ColorSpace', 'srgb')));
nColors = size(RGB_before, 1);
grid_w = 10;
grid_h = nColors / grid_w;
patchSize = 40;

img_before = ones(grid_h * patchSize, grid_w * patchSize, 3);
img_after  = ones(grid_h * patchSize, grid_w * patchSize, 3);

for k = 1:nColors
    row = floor((k-1)/grid_w);
    col = mod((k-1), grid_w);
    r_idx = (row*patchSize+1):((row+1)*patchSize);
    c_idx = (col*patchSize+1):((col+1)*patchSize);
    img_before(r_idx, c_idx, :) = repmat(reshape(RGB_before(k,:),1,1,3), patchSize, patchSize, 1);
    img_after(r_idx, c_idx, :)  = repmat(reshape(RGB_after(k,:),1,1,3), patchSize, patchSize, 1);
end

deltaE = data.deltaE;
dE_grid = zeros(grid_h, grid_w);
for k = 1:nColors
    row = floor((k-1)/grid_w) + 1;
    col = mod((k-1), grid_w) + 1;
    dE_grid(row, col) = deltaE(k);
end

figure('Position', [100 100 1800 600]);
tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
imshow(img_before);
title('Before RGB', 'FontWeight', 'bold', 'FontSize', 15);
hold on;
for k = 1:nColors
    row = floor((k-1)/grid_w);
    col = mod((k-1), grid_w);
    xpos = col*patchSize + patchSize/2;
    ypos = row*patchSize + patchSize/2;
    text(xpos, ypos, num2str(k), 'Color', 'w', 'FontSize', 12, ...
        'HorizontalAlignment','center', 'FontWeight','bold');
end
hold off;

nexttile;
imshow(img_after);
title('After RGB', 'FontWeight', 'bold', 'FontSize', 15);
hold on;
for k = 1:nColors
    row = floor((k-1)/grid_w);
    col = mod((k-1), grid_w);
    xpos = col*patchSize + patchSize/2;
    ypos = row*patchSize + patchSize/2;
    text(xpos, ypos, num2str(k), 'Color', 'w', 'FontSize', 12, ...
        'HorizontalAlignment','center', 'FontWeight','bold');
end
hold off;

nexttile;
imagesc(dE_grid);
axis image off; colormap(jet(255)); colorbar; clim([0 10]);
title('\DeltaE_{2000}', 'FontWeight', 'bold', 'FontSize', 15);
for k = 1:nColors
    row = floor((k-1)/grid_w) + 1;
    col = mod((k-1), grid_w) + 1;
    text(col, row, num2str(k), ...
        'HorizontalAlignment','center', 'Color','w', 'FontSize',12, 'FontWeight','bold');
end

sgtitle('Loaded Palette Patch Order Check', 'FontWeight', 'bold', 'FontSize', 18);
