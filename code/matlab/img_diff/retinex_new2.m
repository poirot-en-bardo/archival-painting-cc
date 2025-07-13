clear; close all; clc;

%% --- Paths & Load ---
path_before = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat';
path_after  = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat';
path_film   = '/home/oem/eliza/data/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat';

b = load(path_before);
a = load(path_after);
f = load(path_film);

Lab_bef    = b.Lab_img;      % before (ground truth)
Lab_aft    = a.Lab_img;      % after
film_Lab   = f.Lab_img;      % film, direct Lab
film_RGB   = f.RGB_lin_img;  % film, linear RGB

%% --- Compute three ΔE maps ---
% 1) direct Film→After in Lab
dE_film_direct = compute_deltaE2000(film_Lab, Lab_aft);

% 2) Retinex Film→After
sigmas = [2, 10, 50];
refl_film = apply_multiscale_retinex(film_RGB, sigmas);
refl_aft  = apply_multiscale_retinex(a.RGB_lin_img, sigmas);
dE_ret    = compute_reflectance_deltaE_direct(refl_film, refl_aft, false);

% 3) Ground-truth Before→After in Lab
dE_gt     = compute_deltaE2000(Lab_bef, Lab_aft);

%% --- Find percentile corresponding to ΔE ≥ 6 on direct Film→After ---
fixed_DE = 6;
pct = 100 * (1 - nnz(dE_film_direct >= fixed_DE) / numel(dE_film_direct));

% Now threshold **all** maps at that percentile
thr_direct = prctile(dE_film_direct(:), pct);
thr_ret    = prctile(dE_ret(:),           pct);
thr_gt     = prctile(dE_gt(:),            pct);  % optional, for consistent display

mask_film_direct = dE_film_direct >= thr_direct;
mask_retinex     = dE_ret           >= thr_ret;
mask_gt          = dE_gt            >= fixed_DE;  % keep GT at ΔE≥6

%% --- Display ---
film_srgb  = xyz2rgb(f.RGB_lin_img,'WhitePoint','d50');
after_srgb = xyz2rgb(a.RGB_lin_img,'WhitePoint','d50');

figure('Name','sRGB Comparison','Units','normalized','Position',[.2 .6 .4 .3]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile; imshow(film_srgb);  title('Film');
nexttile; imshow(after_srgb); title('After');

figure('Name','Binary Masks','Units','normalized','Position',[.1 .1 .8 .4]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile; imshow(mask_film_direct); title(sprintf('Direct ΔE ≥ %.2f', thr_direct));
nexttile; imshow(mask_retinex);     title(sprintf('Retinex ΔE ≥ %.2f', thr_ret));
nexttile; imshow(mask_gt);          title('GT: ΔE_{bef→aft} ≥ 6');

%% --- Helper Functions ---

function dE = compute_deltaE2000(L1,L2)
    sz = size(L1);
    A  = reshape(L1,[],3);
    B  = reshape(L2,[],3);
    d  = deltaE2000(A,B);
    dE = reshape(d, sz(1), sz(2));
end

function R = apply_multiscale_retinex(img, sigmas)
    img = im2double(img); img(img<=0)=eps;
    R = zeros(size(img));
    for s = sigmas
        L = imgaussfilt(rgb2gray(img), s); L(L<=0)=eps;
        for c=1:3
            R(:,:,c) = R(:,:,c) + log(img(:,:,c)) - log(L);
        end
    end
    R = R/numel(sigmas);
    R = R - min(R(:)); R = R / max(R(:));
end

function dE = compute_reflectance_deltaE_direct(r1,r2,isGamma)
    [h,w,~] = size(r1);
    X1 = prophoto2xyz(r1,isGamma);
    X2 = prophoto2xyz(r2,isGamma);
    L1 = xyz2lab_custom(reshape(X1,[],3));
    L2 = xyz2lab_custom(reshape(X2,[],3));
    d  = deltaE2000(L1,L2);
    dE = reshape(d,h,w);
end
