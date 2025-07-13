clear; close all;
addpath('SCIELAB-1996');  % Ensure scielab.m is available

% Set paths
% path_before = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat';
% path_after  = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat';
path_before = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat';
path_after  = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat';

% path_film = '/home/oem/eliza/data/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat';
% path_film = '/home/oem/eliza/data/xyz_lab_rgb/film/cactus_led_fuji_exp0.mat';
% path_film = '/home/oem/eliza/data/xyz_lab_rgb/film/cactus_led_fuji_overexp.mat';
path_film = '/home/oem/eliza/data/xyz_lab_rgb/film/yoda_halogen_fuji_exp0.mat';
% path_film = '/home/oem/eliza/data/xyz_lab_rgb/film/yoda_halogen_fuji_overexp.mat';
% path_film = '/home/oem/eliza/data/xyz_lab_rgb/film/yoda_led_kodak_exp0.mat';

painting_before = load(path_before);
painting_after  = load(path_after);
film_data       = load(path_film);

% Extract original XYZ and gamma-encoded RGB
XYZ_after = painting_after.XYZ_img;
XYZ_film  = film_data.XYZ_img;
RGB_after = painting_after.RGB_img;
RGB_film  = film_data.RGB_img;

% --- Step 1: Retinex Preprocessing on RGB ---
sigma = 5;
refl_after = apply_retinex(RGB_after, sigma);
refl_film  = apply_retinex(RGB_film, sigma);

% --- Step 2: Convert Retinex RGB → XYZ ---
XYZ_refl_after = prophoto2xyz(refl_after, true);  % true = gamma-encoded
XYZ_refl_film  = prophoto2xyz(refl_film,  true);

% --- Step 3: Compute sCIELAB on corrected XYZs ---
dE_map_scielab = compute_scielab_deltaE_from_XYZ(XYZ_refl_film, XYZ_refl_after);

% --- ΔE2000: Ground Truth (HSI before vs after) ---
Lab_before = painting_before.Lab_img;
Lab_after  = painting_after.Lab_img;
dE_map_direct = reshape(deltaE2000(reshape(Lab_before, [], 3), reshape(Lab_after, [], 3)), size(Lab_after,1), size(Lab_after,2));

% --- ΔE2000: Film vs HSI after ---
Lab_film = film_data.Lab_img;
dE_map_film_direct = reshape(deltaE2000(reshape(Lab_film, [], 3), reshape(Lab_after, [], 3)), size(Lab_after,1), size(Lab_after,2));

% --- Thresholding ---
use_fixed_threshold = false;
if use_fixed_threshold
    deltaE_threshold = 6;
    mask_direct       = dE_map_direct       >= deltaE_threshold;
    mask_film_direct  = dE_map_film_direct  >= deltaE_threshold;
    mask_scielab      = dE_map_scielab      >= deltaE_threshold;
else
    p = 70;
    mask_direct       = dE_map_direct      >= prctile(dE_map_direct(:),      p);
    mask_film_direct  = dE_map_film_direct >= prctile(dE_map_film_direct(:), p);
    mask_scielab      = dE_map_scielab     >= prctile(dE_map_scielab(:),     p);
end

% --- Visualization: sRGB images ---
film_srgb  = xyz2rgb(XYZ_film ./ 100, 'WhitePoint', 'd50');
after_srgb = xyz2rgb(XYZ_after ./ 100, 'WhitePoint', 'd50');

figure('Name', 'sRGB Comparison', 'Units', 'normalized', 'Position', [0.2 0.6 0.4 0.5]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
fontSize = 18;
nexttile; imshow(film_srgb);  title('Film (sRGB)', 'FontSize', fontSize); axis off;
nexttile; imshow(after_srgb); title('HSI After (sRGB)', 'FontSize', fontSize); axis off;

% --- Visualization: Thresholded ΔE Maps ---
figure('Name', 'Thresholded ΔE Maps', 'Units', 'normalized', 'Position', [0.05 0.1 0.75 0.6]);
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile; imshow(mask_film_direct); title({'Film vs. HSI after', 'Lab ΔE2000'}, 'FontSize', fontSize); axis off;
nexttile; imshow(mask_scielab);     title({'Film vs. HSI after', 'sCIELAB + Retinex'}, 'FontSize', fontSize); axis off;
nexttile; imshow(mask_direct);      title({'HSI before vs. after', 'Ground Truth'}, 'FontSize', fontSize); axis off;

%% --- Retinex Preprocessing ---
function reflectance = apply_retinex(img_rgb, sigma)
    img_rgb = im2double(img_rgb);
    img_rgb(img_rgb <= 0) = eps;
    illum = imgaussfilt(rgb2gray(img_rgb), sigma);
    illum(illum <= 0) = eps;
    reflectance = zeros(size(img_rgb));
    for c = 1:3
        reflectance(:,:,c) = log(img_rgb(:,:,c)) - log(illum);
    end
    reflectance = reflectance - min(reflectance(:));
    reflectance = reflectance / max(reflectance(:));
end

%% --- sCIELAB in XYZ ---
function dE_map = compute_scielab_deltaE_from_XYZ(XYZ1, XYZ2)
    whitepoint = [0.9642 1.0000 0.8251];  % D50
    ppd = 10;  % pixels per degree
    dE_map = scielab(ppd, XYZ1, XYZ2, whitepoint, 'XYZ');
end
