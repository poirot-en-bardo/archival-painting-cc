clear; close all; clc;

%% --- User‐set threshold parameters ---
% Option A: absolute Mahalanobis-distance cutoff
thr_manual = 2.3;        % e.g. keep pixels with dM ≥ thr_manual

% Option B: percentile threshold
p_manual   = 60;         % e.g. keep the top p_manual% of dM values

%% --- Paths & Load Data ---
path_before = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat';
path_after  = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat';
path_film   = '/home/oem/eliza/data/xyz_lab_rgb/film/yoda_halogen_fuji_exp0.mat';
path_before = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat';
path_after  = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat';
path_film   = '/home/oem/eliza/data/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat';

d_before = load(path_before);
d_after  = load(path_after);
d_film   = load(path_film);

Lab_bef  = d_before.Lab_img;    % M×N×3 before
Lab_aft  = d_after.Lab_img;     %      after
Lab_film = d_film.Lab_img;      %      film

[M, N, ~] = size(Lab_film);

%% --- Ground Truth Mask (ΔE_before→after ≥ 6) ---
dE_gt   = compute_deltaE2000(Lab_bef, Lab_aft);
mask_gt = dE_gt >= 6;

%% --- Compute Fast Mahalanobis‐CVA Distances ---

% 1) Build per-pixel Lab-difference matrix D of size (M*N)×3
D = reshape(Lab_aft - Lab_film, [], 3);    % double

% 2) Sample covariance (fast)
C = cov(D);

% 3) Cholesky factorization: C = L * L'
Lchol = chol(C + 1e-6*eye(3), 'lower');

% 4) Chunked solve for Mahalanobis distances
nPix     = size(D,1);
Dsingle  = single(D);
dMvec    = zeros(nPix,1,'single');
batch    = 200e3;   % ~200k pixels per chunk

for idx0 = 1 : batch : nPix
    idx  = idx0 : min(idx0 + batch - 1, nPix);
    Di   = Dsingle(idx,:);          % b×3
    Yi   = Lchol \ Di';             % 3×b
    dMvec(idx) = sqrt( sum(Yi.^2,1) )';  % b×1
end

dM = reshape(double(dMvec), M, N);  % back to M×N

%% --- Thresholding Options ---

% Absolute threshold mask
mask_abs = dM >= thr_manual;

% Percentile threshold mask
thrP      = prctile(dM(:), p_manual);
mask_pct  = dM >= thrP;

%% --- Display GT, Absolute‐Thresh, Percentile‐Thresh Side by Side ---

figure('Name','Mahalanobis‐CVA Threshold Comparison','Units','normalized','Position',[0.1 0.2 0.8 0.5]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

nexttile;
imshow(mask_gt);
title('GT: ΔE_{bef→aft} ≥ 6');

nexttile;
imshow(mask_abs);
title(sprintf('Mahalanobis ≥ %.2f', thr_manual));

nexttile;
imshow(mask_pct);
title(sprintf('Mahalanobis ≥ %d%% pct', p_manual));

%% --- Helper: Compute ΔE2000 map ---
function dE_map = compute_deltaE2000(L1, L2)
    % L1, L2: M×N×3 Lab images
    sz = size(L1);
    A  = reshape(L1, [], 3);
    B  = reshape(L2, [], 3);
    dE = deltaE2000(A, B);
    dE_map = reshape(dE, sz(1), sz(2));
end
