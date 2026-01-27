%% MTT-based change detection in log-RGB
% Requires: both images loaded as linear RGB in range [0,1] (RGB_lin_img)
% film_data.RGB_lin_img  and painting_after.RGB_lin_img are expected
clear; close all;
% --- Paths (edit to your files) ---
%% --- Paths ---
path_before = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat'; 
path_after = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat'; 
% path_before = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat'; 
% path_after = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat';



% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat';
path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_exp0.mat'; 
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_underexp.mat';
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/yoda_halogen_fuji_exp0.mat';
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/yoda_halogen_fuji_overexp.mat';
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/yoda_led_kodak_exp0.mat';


% oRGB-based change detection (film vs after)
% - Requires: film_data.RGB_lin_img and painting_after.RGB_lin_img (linear RGB in 0..1)
% - Uses oRGB transform (Bratkova / Boulos / Shirley) then computes:
%     * Euclidean distance in oRGB chrominance (rg, yb)
%     * Multi-scale DoG on the oRGB chrominance channels (per-scale Euclidean, mean across scales)
% - Adjust parameters (K, DoG scales, thresholds) as needed.




%% --- Load data ---
painting_before = load(path_before);
painting_after  = load(path_after);
film   = load(path_film);

% Linear RGB in [0,1]
RGB_before = painting_before.RGB_lin_img;
RGB_after  = painting_after.RGB_lin_img;
RGB_film   = film.RGB_lin_img;
RGB_before = painting_before.RGB_img;
RGB_after  = painting_after.RGB_img;
RGB_film   = film.RGB_img;
% 
% RGB_before = xyz2rgb(before.XYZ_img, "ColorSpace","linear-rgb","WhitePoint","d65", "OutputType","double");
% RGB_after  = xyz2rgb(after.XYZ_img, "ColorSpace","linear-rgb","WhitePoint","d65", "OutputType","double");
% RGB_film   = xyz2rgb(film.XYZ_img, "ColorSpace","linear-rgb","WhitePoint","d65", "OutputType","double");

[h,w,~] = size(RGB_film);
min_dim = min(h,w);

%% --- Parameters ---
eps_val = 1e-9;
DoG_small = 1;
sigma_large_list = [min_dim/10, min_dim/5, min_dim/3];
% sigma_large_list = [min_dim/5, min_dim/3, min_dim/2];
nScales = numel(sigma_large_list);
DoG_at = @(X,sL) imgaussfilt(X, DoG_small) - imgaussfilt(X, sL);

normalizeLuma = true;
clipChromRange = true;

%% --- Helper function ---
function theta_o = piecewise_theta_o(theta)
    theta2 = mod(theta + 2*pi, 2*pi);
    theta_o2 = zeros(size(theta2));
    half_mask = theta2 <= pi;
    t_half = theta2(half_mask);
    theta_o2(half_mask) = (3/2) * t_half;
    idx_p2 = half_mask & (theta2 > pi/3);
    t2 = theta2(idx_p2);
    theta_o2(idx_p2) = pi/2 + (3/4)*(t2 - pi/3);
    idx_reflect = theta2 > pi;
    t_ref = 2*pi - theta2(idx_reflect);
    theta_ref_out = zeros(size(t_ref));
    mask_r1 = t_ref <= (pi/3);
    theta_ref_out(mask_r1) = (3/2) * t_ref(mask_r1);
    theta_ref_out(~mask_r1) = pi/2 + (3/4)*(t_ref(~mask_r1)-pi/3);
    theta_o2(idx_reflect) = -theta_ref_out;
    theta_o = theta_o2;
    theta_o(theta_o>pi) = theta_o(theta_o>pi) - 2*pi;
end

%% --- Function to compute chroma channels ---
function [C1,C2] = compute_oRGB_chroma(R,G,B,normalizeLuma,clipChromRange)
    eps_val = 1e-9;
    % Luma
    L = 0.299*R + 0.587*G + 0.114*B;
    yb = 0.5*(R+G)-B;
    rg = 0.866*(R-G);
    
    if normalizeLuma
        Lmin = min(L(:));
        Lmax = max(L(:));
        L = (L-Lmin)/(Lmax-Lmin+eps_val);
    end
    
    % Polar coordinates
    theta = atan2(rg,yb);
    theta_o = piecewise_theta_o(theta);
    delta = theta_o - theta;
    
    % Rotate chroma vectors
    C1 = cos(delta).*yb - sin(delta).*rg;
    C2 = sin(delta).*yb + cos(delta).*rg;
    
    % Clip range
    if clipChromRange
        maxMag = max([max(abs(C1(:))), max(abs(C2(:)))]);
        if maxMag>0
            C1 = C1 / (maxMag+eps_val);
            C2 = C2 / (maxMag+eps_val);
        end
    end
end

%% --- Compute change maps function ---
function [delta_chroma,DoG_fused] = change_maps(C1A,C2A,C1B,C2B,sigma_large_list,DoG_at)
    delta_chroma = sqrt((C1A-C1B).^2 + (C2A-C2B).^2);
    delta_chroma = mat2gray(delta_chroma);
    
    nScales = numel(sigma_large_list);
    [h,w] = size(C1A);
    DoG_maps = zeros(h,w,nScales);
    for si=1:nScales
        sL = sigma_large_list(si);
        d1 = DoG_at(C1A,sL) - DoG_at(C1B,sL);
        d2 = DoG_at(C2A,sL) - DoG_at(C2B,sL);
        DoG_maps(:,:,si) = sqrt(d1.^2 + d2.^2);
    end
    DoG_fused = mat2gray(mean(DoG_maps,3));
end

%% --- Compute chroma channels ---
[C1_F,C2_F] = compute_oRGB_chroma(RGB_film(:,:,1),RGB_film(:,:,2),RGB_film(:,:,3),normalizeLuma,clipChromRange);
[C1_A,C2_A] = compute_oRGB_chroma(RGB_after(:,:,1),RGB_after(:,:,2),RGB_after(:,:,3),normalizeLuma,clipChromRange);
[C1_B,C2_B] = compute_oRGB_chroma(RGB_before(:,:,1),RGB_before(:,:,2),RGB_before(:,:,3),normalizeLuma,clipChromRange);

%% --- Compute change maps ---
[delta_chroma_FA,DoG_FA] = change_maps(C1_F,C2_F,C1_A,C2_A,sigma_large_list,DoG_at);
[delta_chroma_BA,DoG_BA] = change_maps(C1_B,C2_B,C1_A,C2_A,sigma_large_list,DoG_at);

%% --- Visualize ---
figure('Name','Change Maps: Film vs After','Position',[100 100 1200 500]);
subplot(1,2,1); imagesc(delta_chroma_FA); axis image off; colorbar; colormap(turbo); title('Δ chroma: Film vs After');
subplot(1,2,2); imagesc(DoG_FA); axis image off; colorbar; colormap(turbo); title('DoG fused: Film vs After');

figure('Name','Change Maps: Before vs After','Position',[100 100 1200 500]);
subplot(1,2,1); imagesc(delta_chroma_BA); axis image off; colorbar; colormap(turbo); title('Δ chroma: Before vs After');
subplot(1,2,2); imagesc(DoG_BA); axis image off; colorbar; colormap(turbo); title('DoG fused: Before vs After');


%%
%% --- Binary masks based on DoG threshold ---

percentile_val = 80; 

% Before vs After
thresh_BA = prctile(DoG_BA(:), percentile_val);
mask_BA = DoG_BA > thresh_BA;

figure('Name','Binary Mask: Before vs After','Position',[100 100 600 500]);
imagesc(mask_BA); axis image off; colormap(gray); colorbar;
title(['Binary mask: Before vs After (>', num2str(percentile_val), 'th percentile)']);

% Film vs After
thresh_FA = prctile(DoG_FA(:), percentile_val);
mask_FA = DoG_FA > thresh_FA;

figure('Name','Binary Mask: Film vs After','Position',[100 100 600 500]);
imagesc(mask_FA); axis image off; colormap(gray); colorbar;
title(['Binary mask: Film vs After (>', num2str(percentile_val), 'th percentile)']);











%%
%% --------------------------------------------------------------
%  Compute ΔE2000 map (Before vs After) and threshold at 6
% --------------------------------------------------------------

LAB_before = painting_before.Lab_img;
LAB_after  = painting_after.Lab_img;

[h, w, ~] = size(LAB_before);
DE_map = zeros(h, w);

for i = 1:h
    for j = 1:w
        DE_map(i,j) = deltaE2000( squeeze(LAB_before(i,j,:))', ...
                                  squeeze(LAB_after(i,j,:))', ...
                                  [1 1 1] );
    end
end

% Ground-truth mask from ΔE >= 6
GT_mask = DE_map >= 6;

figure; imagesc(GT_mask); axis image off; colormap(gray);
title('Ground-truth ΔE2000 ≥ 6');


%% --------------------------------------------------------------
%  Find the percentile that best matches the ΔE GT mask
%  using DoG-based change map (e.g. DoG_FA)
% --------------------------------------------------------------

P_range = 50:99;     % search percentiles
F1_list = zeros(size(P_range));
Acc_list = zeros(size(P_range));
Prec_list = zeros(size(P_range));
Rec_list = zeros(size(P_range));
IoU_list = zeros(size(P_range));


DoG_map = DoG_FA;    % <<< choose your DoG map here

for idx = 1:numel(P_range)
    p = P_range(idx);

    % Threshold at percentile p
    thresh = prctile(DoG_map(:), p);
    mask_p = DoG_map > thresh;

    % Ensure sizes match
    if ~isequal(size(mask_p), size(GT_mask))
        error('ERROR: mask and GT sizes do not match.');
    end

    % Compute confusion matrix components
    TP = sum(mask_p(:) & GT_mask(:));
    FP = sum(mask_p(:) & ~GT_mask(:));
    FN = sum(~mask_p(:) & GT_mask(:));
    TN = sum(~mask_p(:) & ~GT_mask(:));

    % Metrics
    Prec = TP / (TP + FP + eps);
    Rec  = TP / (TP + FN + eps);
    F1   = 2*(Prec*Rec)/(Prec+Rec+eps);
    Acc  = (TP + TN) / (TP + TN + FP + FN);

    Prec_list(idx) = Prec;
    Rec_list(idx)  = Rec;
    F1_list(idx)   = F1;
    Acc_list(idx)  = Acc;
    IoU_list(idx) = TP / (TP + FP + FN + eps); % Jaccard Index

end

% Find best percentile (by F1 score)
[bestF1, bestIdx] = max(F1_list);
bestPercentile = P_range(bestIdx);
bestAcc = Acc_list(bestIdx);

fprintf('\n==============================\n');
fprintf(' Optimal percentile = %d\n', bestPercentile);
fprintf(' Best F1 = %.4f\n', bestF1);
fprintf(' Accuracy = %.4f\n', bestAcc);
fprintf('==============================\n\n');

%% --------------------------------------------------------------
%  Visualize optimal-threshold binary mask
% --------------------------------------------------------------
opt_thresh = prctile(DoG_map(:), bestPercentile);
best_mask = DoG_map > opt_thresh;

figure; 
subplot(1,2,1); imagesc(GT_mask); axis image off; colormap(gray);
title('Ground Truth ΔE ≥ 6');

subplot(1,2,2); imagesc(best_mask); axis image off; colormap(gray);
title(['Best DoG mask (percentile ' num2str(bestPercentile) ')']);
