clear; close all; clc;

%% --- User parameter: top‐percentile for all methods ---
p_thresh = 60;   % keep only the top (100-p_thresh)% strongest changes

%% --- Paths & Load ---
p1 = load('/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat');
p2 = load('/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat');
pf = load('/home/oem/eliza/data/xyz_lab_rgb/film/yoda_halogen_fuji_exp0.mat');
p1 = load('/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat');
p2  = load('/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat');
pf   = load('/home/oem/eliza/data/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat');

Lab_b = p1.Lab_img;   % M×N×3
Lab_a = p2.Lab_img;
Lab_f = pf.Lab_img;
[M,N,~] = size(Lab_f);

%% --- Ground‐truth mask ΔE(before→after) ≥ 6 ---
dE_gt = compute_deltaE2000(Lab_b, Lab_a);
mask_gt = dE_gt >= 6;

%% --- Precompute film→after ΔE and Lab‐difference ---
dE_fa   = compute_deltaE2000(Lab_f, Lab_a);
diffLab = reshape(Lab_a - Lab_f, M*N, 3);

%% --- Method 1: Mahalanobis CVA ---
C     = cov(diffLab);         % 3×3 covariance
iC    = inv(C + 1e-6*eye(3)); % regularize
% Mahalanobis distance per‐pixel:
dMvec = sqrt( sum((diffLab * iC) .* diffLab, 2) );
dM    = reshape(dMvec, M, N);
thr1  = prctile(dM(:), p_thresh);
mask1 = dM >= thr1;

%% --- Method 2: Difference‐of‐Gaussians on ΔE map ---
small_sigma = 2;
large_sigma = 120;
g1 = imgaussfilt(dE_fa, small_sigma);
g2 = imgaussfilt(dE_fa, large_sigma);
dog = abs(g1 - g2);
thr2 = prctile(dog(:), p_thresh);
mask2 = dog >= thr2;

%% --- Method 3: Edge‐XOR on L channel ---
L_f = mat2gray(Lab_f(:,:,1));
L_a = mat2gray(Lab_a(:,:,1));
e1  = edge(L_f, 'Canny');
e2  = edge(L_a, 'Canny');
mask3 = xor(e1, e2);

%% --- Display all 4 masks ---
figure('Units','normalized','Position',[.1 .1 .8 .8]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile; imshow(mask_gt); title('GT: ΔE_{bef→aft}≥6');
nexttile; imshow(mask1);   title(sprintf('Mahalanobis CVA ≥ %d%%', p_thresh));
nexttile; imshow(mask2);   title(sprintf('DoG‐ΔE ≥ %d%%', p_thresh));
nexttile; imshow(mask3);   title('Edge‐XOR (Canny)');

%% --- Helper: ΔE2000 map ---
function M = compute_deltaE2000(L1, L2)
    sz = size(L1);
    A  = reshape(L1, [], 3);
    B  = reshape(L2, [], 3);
    d  = deltaE2000(A, B);
    M  = reshape(d, sz(1), sz(2));
end
