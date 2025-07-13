clear; close all; clc;

%% --- Paths & Load ---
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
dE_gt    = compute_deltaE2000(Lab_bef, Lab_aft);
mask_gt  = dE_gt >= 6;
frac_gt  = nnz(mask_gt) / numel(mask_gt);   % fraction of “changed” pixels

%% --- Improved Mahalanobis‐CVA ---

% 1) Build per-pixel Lab difference vectors for film→after
D = reshape(Lab_aft - Lab_film, [], 3);   % (M*N)×3

% 2) Robust covariance estimate (less sensitive to outliers)
%    Requires Statistics Toolbox.  If unavailable, use cov(D).
try
    [covRob, ~] = robustcov(D);
catch
    covRob = cov(D);
end

% 3) Cholesky factorization (covRob = Lchol * Lchol')
Lchol = chol(covRob + 1e-6*eye(3), 'lower');

% 4) Solve Lchol * Y = D'  →  Y = Lchol \ D'
Y = Lchol \ D';   % 3×(M*N)

% 5) Mahalanobis distance per pixel
dMvec = sqrt(sum(Y.^2, 1))';    % (M*N)×1
dM    = reshape(dMvec, M, N);   % back to M×N

% 6) Automatic threshold: pick the same fraction as GT
thrM  = prctile(dM(:), 100*(1 - frac_gt));
mask_m = dM >= thrM;

%% --- Display Results ---
figure('Name','Mahalanobis‐CVA vs ΔE Ground Truth','Units','normalized','Position',[0.1 0.2 0.6 0.6]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile;
imshow(mask_gt);
title('GT: ΔE_{bef→aft} ≥ 6');

nexttile;
imshow(mask_m);
title('Mahalanobis‐CVA (robust cov, matched frac)');

%% --- Helper: Compute ΔE2000 map ---
function dE_map = compute_deltaE2000(Lab1, Lab2)
    % Lab1, Lab2: M×N×3
    sz = size(Lab1);
    A  = reshape(Lab1, [], 3);
    B  = reshape(Lab2, [], 3);
    dE = deltaE2000(A, B);
    dE_map = reshape(dE, sz(1), sz(2));
end
