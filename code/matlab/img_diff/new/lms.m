%% --- Paths ---
path_before = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat'; 
path_after = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat'; 
path_before = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat'; 
path_after = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat';



% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat';
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_exp0.mat';
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_underexp.mat';
path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/yoda_halogen_fuji_exp0.mat';
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/yoda_halogen_fuji_overexp.mat';
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/yoda_led_kodak_exp0.mat';





%% --- Load data ---
painting_after  = load(path_after);
film_data       = load(path_film);

XYZ_after  = painting_after.XYZ_img;
XYZ_film   = film_data.XYZ_img;

[h,w,~] = size(XYZ_film);
min_dim = min(h,w);

%% --- Parameters ---
sigma_large_list = [min_dim/10, min_dim/5, min_dim/3];
nScales = numel(sigma_large_list);
DoG_small = 1;
DoG_at = @(X, sL) imgaussfilt(X, DoG_small) - imgaussfilt(X, sL);

eps_val = 1e-9; % avoid division by zero

%% --- XYZ → LMS conversion (Stockman-Sharpe 2° cone fundamentals, D65) ---
% Matrix from XYZ to LMS (approximate)
M_XYZ2LMS = [ ...
    0.4002 0.7075 -0.0808; ...
   -0.2263 1.1653  0.0457; ...
    0       0      0.9182];

% Reshape for matrix multiplication
F_film  = reshape(XYZ_film, [], 3)';
F_after = reshape(XYZ_after, [], 3)';

LMS_film  = M_XYZ2LMS * F_film;
LMS_after = M_XYZ2LMS * F_after;

% Reshape back to h x w x 3
LMS_film  = reshape(LMS_film', h, w, 3);
LMS_after = reshape(LMS_after', h, w, 3);

%% --- Compute opponent channels ---
L_film = LMS_film(:,:,1); M_film = LMS_film(:,:,2); S_film = LMS_film(:,:,3);
L_after = LMS_after(:,:,1); M_after = LMS_after(:,:,2); S_after = LMS_after(:,:,3);

% Red-Green (O1) and Blue-Yellow (O2)
O1_film  = L_film - M_film;
O2_film  = L_film + M_film - S_film;

O1_after = L_after - M_after;
O2_after = L_after + M_after - S_after;

% Optional normalization by total intensity
tot_film  = L_film + M_film + S_film + eps_val;
tot_after = L_after + M_after + S_after + eps_val;

O1_film  = O1_film ./ tot_film;
O2_film  = O2_film ./ tot_film;
O1_after = O1_after ./ tot_after;
O2_after = O2_after ./ tot_after;

%% --- Compute ΔO map ---
deltaO = sqrt((O1_film - O1_after).^2 + (O2_film - O2_after).^2);
deltaO_norm = mat2gray(deltaO);

%% --- Compute DoG map on opponent channels ---
DoG_maps = zeros(h,w,nScales);
for si = 1:nScales
    sL = sigma_large_list(si);
    dO1 = DoG_at(O1_film, sL) - DoG_at(O1_after, sL);
    dO2 = DoG_at(O2_film, sL) - DoG_at(O2_after, sL);
    DoG_maps(:,:,si) = sqrt(dO1.^2 + dO2.^2);
end
DoG_opponent = mat2gray(sum(DoG_maps,3));

%% --- Visualization ---
figure('Name','ΔO1O2 Map (LMS Opponent)','Position',[50 50 800 400]);
imagesc(deltaO_norm); axis image off; colorbar; colormap(jet);
title('ΔO1O2 (LMS Opponent Channels)');

figure('Name','DoG Opponent Map (LMS)','Position',[50 50 800 400]);
imagesc(DoG_opponent); axis image off; colorbar; colormap(turbo);
title('DoG Opponent Channels');

%% --- Thresholded binary maps ---
threshold_delta = 0.05;
threshold_DoG   = 0.05;

figure('Name','Binary ΔO Map','Position',[50 50 800 400]);
imshow(deltaO_norm >= threshold_delta); title('ΔO1O2 ≥ threshold');

figure('Name','Binary DoG Opponent Map','Position',[50 50 800 400]);
imshow(DoG_opponent >= threshold_DoG); title('DoG Opponent ≥ threshold');
