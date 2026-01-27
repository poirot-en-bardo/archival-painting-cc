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

%% --- Compute opponent color ratios ---
% First opponent channels (O1, O2)
O1_film  = (XYZ_film(:,:,1) - XYZ_film(:,:,3)) ./ (XYZ_film(:,:,1) + XYZ_film(:,:,3) + XYZ_film(:,:,2) + eps_val);
O2_film  = (2*XYZ_film(:,:,2) - XYZ_film(:,:,1) - XYZ_film(:,:,3)) ./ (XYZ_film(:,:,1) + XYZ_film(:,:,3) + XYZ_film(:,:,2) + eps_val);

O1_after = (XYZ_after(:,:,1) - XYZ_after(:,:,3)) ./ (XYZ_after(:,:,1) + XYZ_after(:,:,3) + XYZ_after(:,:,2) + eps_val);
O2_after = (2*XYZ_after(:,:,2) - XYZ_after(:,:,1) - XYZ_after(:,:,3)) ./ (XYZ_after(:,:,1) + XYZ_after(:,:,3) + XYZ_after(:,:,2) + eps_val);

%% --- Compute ΔO1O2 map ---
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
figure('Name','ΔO1O2 Map','Position',[50 50 800 400]);
imagesc(deltaO_norm); axis image off; colorbar; colormap(turbo);
title('ΔO1O2 (Opponent Ratios)');

figure('Name','DoG Opponent Map','Position',[50 50 800 400]);
imagesc(DoG_opponent); axis image off; colorbar; colormap(turbo);
title('DoG Opponent Ratios');

%% --- Thresholded binary maps ---
threshold_delta = 0.05;
threshold_DoG   = 0.05;

figure('Name','Binary ΔO1O2 Map','Position',[50 50 800 400]);
imshow(deltaO_norm >= threshold_delta); title('ΔO1O2 ≥ threshold');

figure('Name','Binary DoG Opponent Map','Position',[50 50 800 400]);
imshow(DoG_opponent >= threshold_DoG); title('DoG Opponent ≥ threshold');
