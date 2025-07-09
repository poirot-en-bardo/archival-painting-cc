clear; close all;

%% --- INPUT ---
hdrFile = '/home/oem/eliza/data/reflectance/after/cactus_reflectance_after_full.hdr';
outputFolder = '/home/oem/eliza/data/xyz_lab_rgb/reference';

if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

% Load hyperspectral cube
hcube = hypercube(hdrFile);
cube = hcube.DataCube;       % [H x W x B]
wl   = hcube.Wavelength;

% Clip reflectance > 1
cube(cube > 1) = 1;

[H, W, B] = size(cube);

%% --- Load CMFs and D50 Illuminant ---
cmf = importdata('../../../data/CIE2degCMFs_full.csv');
d50 = importdata('../../../data/CIE_D50.txt');
cmf_interp = interp1(cmf(:,1), cmf(:,2:4), wl, 'linear', 'extrap');
ill_interp = interp1(d50(:,1), d50(:,2), wl, 'linear', 'extrap');
S = cmf_interp .* ill_interp;
k = 100 / sum(S(:,2));  % normalization factor

%% --- Generate sRGB preview ---
refl = reshape(cube, [], B);              % [N x B]
XYZ = k * (refl * S);                     % [N x 3]
XYZ_img = reshape(XYZ, H, W, 3);
sRGB = xyz2rgb(XYZ_img ./ 100, 'WhitePoint', 'd50');
sRGB = min(max(sRGB, 0), 1);  % clip

figure('Name', 'sRGB Preview');
imshow(sRGB);
title('sRGB Preview for Patch Selection');

%% --- Draw 24 patches ---
numPatches = 24;
positions = zeros(numPatches, 4);
patchRefl = zeros(numPatches, B);

disp('Draw 24 rectangles around patches...');

for i = 1:numPatches
    h = drawrectangle('Label', sprintf('%d', i), 'Color', 'r');
    wait(h);
    pos = round(h.Position);
    positions(i,:) = pos;

    x1 = max(1, pos(1)); y1 = max(1, pos(2));
    x2 = min(W, x1 + pos(3) - 1); y2 = min(H, y1 + pos(4) - 1);
    patchCube = reshape(cube(y1:y2, x1:x2, :), [], B);
    patchRefl(i,:) = mean(patchCube, 1);  % Mean reflectance per patch
end

%% --- Colorimetric conversion ---
patchXYZ = k * (patchRefl * S);                       % [24 x 3]
patchLab = xyz2lab_custom(patchXYZ);                  % [24 x 3]
patchRGB = xyz2prophoto(patchXYZ ./ 100, true);       % gamma-encoded
patchRGB_lin = xyz2prophoto(patchXYZ ./ 100, false);  % linear

%% --- Save results ---
[~, name, ~] = fileparts(hdrFile);
save(fullfile(outputFolder, [name '_colorchecker_reference.mat']), ...
    'patchRefl', 'patchXYZ', 'patchLab', 'patchRGB', 'patchRGB_lin', 'positions');

fprintf('Saved patch-based color data to:\n%s\n', ...
    fullfile(outputFolder, [name '_colorchecker_reference.mat']));

%% --- Display 6x4 sRGB grid of selected patches ---


patchRGB_srgb = xyz2rgb(patchXYZ ./100, 'ColorSpace', 'srgb', 'WhitePoint', 'd50');
patchRGB_srgb = min(max(patchRGB_srgb, 0), 1);  % Clamp to [0,1]
% patchRGB_srgb = patchRGB;
% Parameters
rows = 4;
cols = 6;
patchSize = 80;  % Size of each square patch (pixels)
numPatches = 24;
% Create grid image
gridImg = zeros(patchSize * rows, patchSize * cols, 3);

for i = 1:numPatches
    r = floor((i-1)/cols);  % row index
    c = mod((i-1), cols);   % col index

    color = reshape(patchRGB_srgb(i,:), 1, 1, 3);
    block = repmat(color, patchSize, patchSize);

    rowStart = r * patchSize + 1;
    colStart = c * patchSize + 1;
    gridImg(rowStart:rowStart+patchSize-1, colStart:colStart+patchSize-1, :) = block;
end

gridImg = min(max(gridImg, 0), 1);  % Just in case

% Show and label
figure('Name', 'Patch RGB Grid (sRGB)');
imshow(gridImg);
title('Reference Target Patches (sRGB)', 'FontSize', 14);

