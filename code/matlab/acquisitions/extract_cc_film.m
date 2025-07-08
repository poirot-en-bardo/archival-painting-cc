clear; close all;

% --- Input MAT file with wb_cube and wl ---
inputFile = '/home/oem/eliza/data/film_scans/white_balanced/yoda_led_kodak_exp0_balanced.mat';
S = load(inputFile);
cube = S.wb_cube;
wl   = S.wl;

% --- Load spectral data ---
D50_path      = '../../../data/CIE_D50.txt';
CMF_path_full = '../../../data/CIE2degCMFs_full.csv';
led_path      = '../../../data/film/CVCL10bands.txt';

D50      = importdata(D50_path);             % [λ, SPD]
CMFs     = importdata(CMF_path_full);        % [λ, x y z]
LEDset   = readmatrix(led_path, 'Delimiter','\t');
led_wl   = LEDset(:,1);
led_spds = LEDset(:,2:end);                  % [L × B]
bandN    = size(led_spds,2);

% --- LED band integration ---
d50_spd = interp1(D50(:,1), D50(:,2), led_wl, 'linear','extrap');
cmf_xyz = interp1(CMFs(:,1), CMFs(:,2:4), led_wl, 'linear','extrap');
d50_band = zeros(bandN,1);
cmf_band = zeros(bandN,3);
for k = 1:bandN
    norm_spd = led_spds(:,k) / sum(led_spds(:,k));
    d50_band(k)   = sum(d50_spd .* norm_spd);
    cmf_band(k,:) = sum(cmf_xyz .* norm_spd,1);
end
data_k_norm = sum(d50_band .* cmf_band(:,2));
k = 100 / data_k_norm;                         % scales so Y = 100

% --- Render to sRGB ---
radCube  = cube .* reshape(d50_band, 1, 1, bandN);
linCube  = reshape(radCube, [], bandN) * cmf_band;
XYZ_img  = reshape(k * linCube, size(cube,1), size(cube,2), 3);
RGB_img  = xyz2rgb(XYZ_img./100, 'ColorSpace','srgb', 'WhitePoint','d50');
RGB_img  = min(max(RGB_img, 0), 1);
%
% --- Show image and select 24 patches
figure('Name','Select 24 Patches'); imshow(RGB_img); title('Select 24 patches (1-by-1)');
%%
numPatches = 24;
patchXYZ = zeros(numPatches, 3);
positions = zeros(numPatches, 4);

for i = 1:numPatches
    h = drawrectangle('Label', sprintf('%d', i), 'Color', 'r');
    wait(h);
    pos = round(h.Position);
    positions(i,:) = pos;

    % Extract mean spectrum from patch
    x1 = max(1, pos(1)); y1 = max(1, pos(2));
    x2 = min(size(cube,2), x1 + pos(3) - 1);
    y2 = min(size(cube,1), y1 + pos(4) - 1);
    patch = cube(y1:y2, x1:x2, :);
    patchAvg = squeeze(mean(reshape(patch, [], bandN), 1));  % [1 x B]

    % Compute XYZ via band-averaged integration
    XYZ = k * (patchAvg .* d50_band') * cmf_band;
    patchXYZ(i,:) = XYZ;
end
%%
% --- Convert to Lab and RGB spaces
patchLab       = xyz2lab(patchXYZ);
patchRGB       = xyz2prophoto(patchXYZ ./100, true);
patchRGB_lin   = xyz2prophoto(patchXYZ ./100, false);  % linear ROMM RGB

% --- Save output
[pathstr, name, ~] = fileparts(inputFile);
outputDir = "/home/oem/eliza/data/film_scans/colorchecker";
outputFile = fullfile(outputDir, [name '_colorchecker.mat']);
save(outputFile, 'patchXYZ', 'patchLab', 'patchRGB', 'patchRGB_lin', 'positions');

fprintf('Saved color checker data to:\n%s\n', outputFile);


%% --- Visualize the 6×4 grid of selected patch RGBs
patchRGB_srgb = xyz2rgb(patchXYZ ./100, 'ColorSpace', 'srgb', 'WhitePoint', 'd50');
patchRGB_srgb = min(max(patchRGB_srgb, 0), 1);  % Clamp to [0,1]
% patchRGB_srgb = patchRGB;
% Parameters
rows = 4;
cols = 6;
patchSize = 50;  % Size of each square patch (pixels)
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
title('Selected Patch Colors (sRGB)', 'FontSize', 14);

