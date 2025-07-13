clear; close all; clc;

%% --- User parameter: CLAHE & Threshold ---
clipLimit   = 0.01;     % between 0 (aggressive) and 1 (none)
numTiles    = [8 8];    % grid of local regions
p_thresh    = 80;       % keep top (100–p_thresh)% of ΔE values

%% --- Paths & Load ---
path_before = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat';
path_after  = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat';
path_film   = '/home/oem/eliza/data/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat';

d_before = load(path_before);
d_after  = load(path_after);
d_film   = load(path_film);

Lab_before = d_before.Lab_img;           % M×N×3
Lab_after  = d_after.Lab_img;            % M×N×3
film_lin   = d_film.RGB_lin_img;         % M×N×3 linear RGB

%% --- Proper CLAHE on L* channel ---
% Convert linear RGB → XYZ → Lab
XYZ       = rgb2xyz(film_lin, 'WhitePoint', 'D50');
Lab_film  = xyz2lab(XYZ,  'WhitePoint', 'D50');  % L* in [0,100]

% Normalize L* to [0,1], apply CLAHE, then restore scale
L  = Lab_film(:,:,1) / 100;
L_eq = adapthisteq(L, ...
        'ClipLimit',     clipLimit, ...
        'NumTiles',      numTiles, ...
        'Distribution',  'uniform');
Lab_film_eq = Lab_film;
Lab_film_eq(:,:,1) = L_eq * 100;

% Convert equalized Lab back to linear RGB for ΔE
XYZ2      = lab2xyz(Lab_film_eq, 'WhitePoint', 'D50');
film_eq_rgb = xyz2rgb(XYZ2,       'WhitePoint', 'D50');

% Convert film_eq_rgb → Lab for ΔE computation
XYZ_eq = rgb2xyz(film_eq_rgb, 'WhitePoint', 'D50');
Lab_eq = xyz2lab(XYZ_eq,       'WhitePoint', 'D50');

%% --- Compute ΔE maps ---
dE_clahe     = compute_deltaE2000(Lab_eq,   Lab_after);
dE_groundtr  = compute_deltaE2000(Lab_before, Lab_after);

%% --- Percentile threshold on CLAHE-based ΔE ---
thr_clahe = prctile(dE_clahe(:), p_thresh);
mask_clahe = dE_clahe >= thr_clahe;

%% --- Ground‐truth mask (ΔE_before→after ≥ 6) ---
mask_gt = dE_groundtr >= 6;

%% --- Cleanup small specks & fill holes ---
minArea = 30;
mask_clahe = bwareaopen(mask_clahe, minArea);
mask_clahe = imfill(mask_clahe, 'holes');
mask_gt     = bwareaopen(mask_gt,     minArea);
mask_gt     = imfill(mask_gt,     'holes');

%% --- Display Results ---
figure('Units','normalized','Position',[0.1 0.1 0.8 0.6]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile;
imshow(mask_clahe);
title(sprintf('CLAHE‐L* ΔE ≥ %d%% pct (thr=%.2f)', p_thresh, thr_clahe));

nexttile;
imshow(mask_gt);
title('GT: ΔE_{before→after} ≥ 6');

%% --- Helper: ΔE2000 computation ---
function dE_map = compute_deltaE2000(Lab1, Lab2)
    sz = size(Lab1);
    A  = reshape(Lab1, [], 3);
    B  = reshape(Lab2, [], 3);
    dE = deltaE2000(A, B);
    dE_map = reshape(dE, sz(1), sz(2));
end
