clear; close all; clc;

%% --- Paths & Load Data ---
path_before = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat';
path_after  = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat';
% path_before = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat';
% path_after  = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat';

path_film = '/home/oem/eliza/data/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat';
% path_film = '/home/oem/eliza/data/xyz_lab_rgb/film/cactus_led_fuji_exp0.mat';
% path_film = '/home/oem/eliza/data/xyz_lab_rgb/film/cactus_led_fuji_underexp.mat';
% path_film = '/home/oem/eliza/data/xyz_lab_rgb/film/yoda_halogen_fuji_exp0.mat';
% path_film = '/home/oem/eliza/data/xyz_lab_rgb/film/yoda_halogen_fuji_overexp.mat';
% path_film = '/home/oem/eliza/data/xyz_lab_rgb/film/yoda_led_kodak_exp0.mat';

d_before = load(path_before);
d_after  = load(path_after);
d_film   = load(path_film);

%% --- Prepare output directory & filename prefix ---
baseDir   = '/home/oem/eliza/masters-thesis/results/plots/dog';
[~, filmName, ~] = fileparts(path_film);
outputDir = fullfile(baseDir, filmName);
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

%% --- Extract RGB images ---
rgb_before = d_before.RGB_img;   % calibrated from MSI scan
rgb_after  = d_after.RGB_img;    % calibrated from HSI scan
rgb_film   = d_film.RGB_img;     % film photograph

%% --- Multiscale High‐Pass Filter in RGB space ---
[M,N,~]      = size(rgb_film);
sigma_min    = 1;
sigma_max    = min(M,N)/2;
sigma_mid    = sqrt(sigma_min * sigma_max);
sigma_scales = [sigma_min, sigma_mid, sigma_max];

hp_film  = multiscale_hp_filter(rgb_film,  sigma_scales);
hp_after = multiscale_hp_filter(rgb_after, sigma_scales);

%% --- Convert HP‐filtered RGB → Lab for ΔE calculation ---
lab_film_hp  = rgb2lab(hp_film,  'WhitePoint','d50');
lab_after_hp = rgb2lab(hp_after, 'WhitePoint','d50');

%% --- Compute ΔE2000 maps ---
% DoG-based change
dE_map_hp     = compute_deltaE2000(lab_film_hp,  lab_after_hp);
% Ground-truth from original Lab data
lab_before = d_before.Lab_img;
lab_after  = d_after.Lab_img;
dE_map_gt    = compute_deltaE2000(lab_before,    lab_after);

%% --- Threshold and create binary masks (ΔE ≥ 6) ---
threshold = 6;
mask_hp   = dE_map_hp  >= threshold;
mask_gt   = dE_map_gt  >= threshold;

%% --- Display results ---
h = figure('Units','normalized','Position',[0.1 0.1 0.8 0.6]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

nexttile;
imshow(rgb_film);
title('Film Photo (RGB)');
axis off;

nexttile;
imshow(mask_hp);
title('DoG‐HP ΔE ≥ 6');
axis off;

nexttile;
imshow(mask_gt);
title('GT ΔE ≥ 6');
axis off;

% Export figure
outFig = fullfile(outputDir, sprintf('%s_DoG_HP_vs_GT.png', filmName));
% exportgraphics(h, outFig, 'Resolution',300);

%% --- Save binary mask for later use ---
maskFile = fullfile(outputDir, sprintf('%s_DoG_HP_Mask.mat', filmName));
mask = logical(mask_hp);
% save(maskFile, 'mask');

%% --- Helper: multiscale high‐pass filter (DoG) ---
function img_hp = multiscale_hp_filter(img_rgb, sigma_scales)
    img_rgb = im2double(img_rgb);
    n = numel(sigma_scales);
    img_hp = zeros(size(img_rgb));
    for c = 1:3
        acc = zeros(size(img_rgb(:,:,1)));
        for sigma = sigma_scales
            blur = imgaussfilt(img_rgb(:,:,c), sigma);
            acc  = acc + (img_rgb(:,:,c) - blur);
        end
        img_hp(:,:,c) = acc / n;
    end
end

%% --- Helper: compute ΔE2000 map ---
function dE_map = compute_deltaE2000(L1, L2)
    sz   = size(L1);
    A    = reshape(L1,[],3);
    B    = reshape(L2,[],3);
    dE   = deltaE2000(A,B);
    dE_map = reshape(dE, sz(1), sz(2));
end
