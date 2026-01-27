clear; close all;

%% --- Paths ---
path_before = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat'; 
path_after = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat'; 
% path_before = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat'; 
% path_after = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat';



% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat';
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_exp0.mat';
path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_underexp.mat';
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/yoda_halogen_fuji_exp0.mat';
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/yoda_halogen_fuji_overexp.mat';
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/yoda_led_kodak_exp0.mat';

%% --- Load data ---
before = load(path_before);
after  = load(path_after);
film   = load(path_film);

Lab_before = before.Lab_img;
Lab_after  = after.Lab_img;
Lab_film   = film.Lab_img;

[h,w,~] = size(Lab_film);
min_dim  = min(h,w);

%% --- Parameters ---
eps_val = 1e-9;
sigma_small = 1;
sigma_large_list = [min_dim/10, min_dim/5, min_dim/3];
nScales = numel(sigma_large_list);

DoG_at = @(X,sL) imgaussfilt(X, sigma_small) - imgaussfilt(X, sL);

%% --- LAB → LCH conversion function ---
function [L, C, h] = Lab_to_LCh(Lab)
    L = Lab(:,:,1);

    a = Lab(:,:,2);
    b = Lab(:,:,3);

    C = sqrt(a.^2 + b.^2);                  % Chroma magnitude
    h = atan2(b, a);                        % Hue angle (radians)
    h = mod(h, 2*pi);                       % 0–2π range
end

%% --- CHANGE MAP FUNCTION ---
function [delta_chroma, DoG_fused] = change_maps_LCh(CA, hA, CB, hB, sigma_large_list, DoG_at)

    % --- chroma diff ---
    deltaC = abs(CA - CB);

    % --- hue diff (circular) ---
    deltaH = abs(atan2(sin(hA-hB), cos(hA-hB)));

    delta_chroma = mat2gray(sqrt(deltaC.^2 + deltaH.^2));

    % --- DoG fused ---
    nScales = numel(sigma_large_list);
    [hgt, wdt] = size(CA);

    DoG_maps = zeros(hgt, wdt, nScales);

    for s = 1:nScales
        sL = sigma_large_list(s);

        d1 = DoG_at(CA, sL) - DoG_at(CB, sL);
        d2 = DoG_at(hA, sL) - DoG_at(hB, sL);

        DoG_maps(:,:,s) = sqrt(d1.^2 + d2.^2);
    end

    DoG_fused = mat2gray(mean(DoG_maps,3));
end

%% --- Convert ALL images to LCh ---
[L_F, C_F, h_F] = Lab_to_LCh(Lab_film);
[L_A, C_A, h_A] = Lab_to_LCh(Lab_after);
[L_B, C_B, h_B] = Lab_to_LCh(Lab_before);

%% --- Compute change maps ---

% Film vs After
[deltaC_FA, DoG_FA] = change_maps_LCh(C_F, h_F, C_A, h_A, sigma_large_list, DoG_at);

% Before vs After
[deltaC_BA, DoG_BA] = change_maps_LCh(C_B, h_B, C_A, h_A, sigma_large_list, DoG_at);

%% --- Visualize Film vs After ---
figure('Name','LCh Change Maps: Film vs After','Position',[100 100 1200 500]);
subplot(1,2,1);
imagesc(deltaC_FA); axis image off; colorbar; colormap(turbo);
title('Δ chroma (LCh): Film vs After');

subplot(1,2,2);
imagesc(DoG_FA); axis image off; colorbar; colormap(turbo);
title('DoG fused (LCh): Film vs After');

%% --- Visualize Before vs After ---
figure('Name','LCh Change Maps: Before vs After','Position',[100 100 1200 500]);
subplot(1,2,1);
imagesc(deltaC_BA); axis image off; colorbar; colormap(turbo);
title('Δ chroma (LCh): Before vs After');

subplot(1,2,2);
imagesc(DoG_BA); axis image off; colorbar; colormap(turbo);
title('DoG fused (LCh): Before vs After');

