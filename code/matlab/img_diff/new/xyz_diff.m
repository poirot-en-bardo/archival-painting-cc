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

close all;

painting_before = load(path_before);
painting_after  = load(path_after);
film_data       = load(path_film);

img_before = painting_before.RGB_img;
img_after   = painting_after.RGB_img;
film_rgb    = film_data.RGB_img;
hyper_rgb = img_after;

%
% ground truth
Lab_before = painting_before.Lab_img;
Lab_after  = painting_after.Lab_img;
Lab_film = film_data.Lab_img;
Lab_hyper = Lab_after;
% dE_map_direct = reshape(deltaE2000(reshape(Lab_before,[],3), reshape(Lab_after,[],3)), size(Lab_before,1), size(Lab_before,2));






%% here
%% --- Parameters ---
scales = [15, 80, 250];  % MSR scales
sigma_small = 1;
sigma_large = 300;
DoG = @(X) imgaussfilt(X, sigma_small) - imgaussfilt(X, sigma_large);

%% --- Load XYZ images ---
XYZ_film  = film_data.XYZ_img;       % HxWx3
XYZ_hyper = painting_after.XYZ_img;  % HxWx3
XYZ_film = prophoto2xyz(film_data.RGB_img, true);
XYZ_hyper = prophoto2xyz(painting_after.RGB_img, true);
[h,w,~] = size(XYZ_film);

%% --- Prepare 3 MSR cases using apply_msr_to_Y ---
film_XYZ_case1  = XYZ_film;                       % no MSR
hyper_XYZ_case1 = XYZ_hyper;

film_XYZ_case2  = apply_msr_to_Y(XYZ_film,  scales); % both MSR
hyper_XYZ_case2 = apply_msr_to_Y(XYZ_hyper, scales);

film_XYZ_case3  = apply_msr_to_Y(XYZ_film, scales);  % film only
hyper_XYZ_case3 = XYZ_hyper;

cases = {
    {'No MSR',        film_XYZ_case1, hyper_XYZ_case1};
    {'Both MSR',      film_XYZ_case2, hyper_XYZ_case2};
    {'Film only MSR', film_XYZ_case3, hyper_XYZ_case3};
};

%% --- Compute ΔX/Z and DoG differences ---
deltaXZ_all = cell(3,1);
DoGXZ_all   = cell(3,1);

for k = 1:3
    XYZF = cases{k}{2};
    XYZH = cases{k}{3};

    % Extract X and Z channels
    XF = XYZF(:,:,1);
    ZF = XYZF(:,:,3);
    XH = XYZH(:,:,1);
    ZH = XYZH(:,:,3);

    % Compute Euclidean difference in XZ plane
    deltaXZ_all{k} = mat2gray(sqrt((XF - XH).^2 + (ZF - ZH).^2));

    % Compute DoG differences in XZ plane
    DoGXZ_all{k} = mat2gray(sqrt((DoG(XF) - DoG(XH)).^2 + (DoG(ZF) - DoG(ZH)).^2));
end

%% --- Visualization ---
figure('Name','ΔX/Z and DoG XZ (MSR cases)','Position',[100 100 1200 500]);
for k = 1:3
    subplot(2,3,k); imagesc(deltaXZ_all{k}); axis image off; colorbar; colormap(jet);
    title(['ΔX/Z – ' cases{k}{1}]);
    subplot(2,3,3+k); imagesc(DoGXZ_all{k}); axis image off; colorbar; colormap(turbo);
    title(['DoG XZ – ' cases{k}{1}]);
end








%% --- Compute full XYZ Euclidean differences and DoG differences ---
deltaXYZ_all = cell(3,1);
DoGXYZ_all   = cell(3,1);
DoG = @(X) imgaussfilt(X, sigma_small) - imgaussfilt(X, sigma_large);

for k = 1:3
    XYZF = cases{k}{2};
    XYZH = cases{k}{3};
    
    % --- Full XYZ Euclidean difference ---
    deltaXYZ_all{k} = mat2gray(sqrt(sum((XYZF - XYZH).^2, 3)));
    
    % --- DoG-based difference in XYZ ---
    DoGXYZ_all{k} = mat2gray(sqrt(sum((DoG(XYZF) - DoG(XYZH)).^2, 3)));
end

%% --- Visualization ---
figure('Name','Full XYZ and DoG XYZ Differences','Position',[100 100 1200 600]);
for k = 1:3
    subplot(2,3,k);
    imagesc(deltaXYZ_all{k});
    axis image off;
    colorbar;
    colormap(turbo);
    title(['ΔXYZ – ' cases{k}{1}]);
    
    subplot(2,3,3+k);
    imagesc(DoGXYZ_all{k});
    axis image off;
    colorbar;
    colormap(turbo);
    title(['DoG XYZ – ' cases{k}{1}]);
end







%%
%% --- Parameters ---
sigma_small = 1;
sigma_large = 300;
DoG = @(X) imgaussfilt(X, sigma_small) - imgaussfilt(X, sigma_large);

%% --- Load XYZ images ---
XYZ_before = painting_before.XYZ_img;  % HxWx3
XYZ_after  = painting_after.XYZ_img;   % HxWx3

%% --- Compute ΔXZ, ΔXYZ and DoG variants ---
deltaXZ  = sqrt((XYZ_before(:,:,1) - XYZ_after(:,:,1)).^2 + ...
                (XYZ_before(:,:,3) - XYZ_after(:,:,3)).^2);

deltaXYZ = sqrt(sum((XYZ_before - XYZ_after).^2, 3));

DoGXZ  = sqrt((DoG(XYZ_before(:,:,1)) - DoG(XYZ_after(:,:,1))).^2 + ...
              (DoG(XYZ_before(:,:,3)) - DoG(XYZ_after(:,:,3))).^2);

DoGXYZ = sqrt(sum((DoG(XYZ_before) - DoG(XYZ_after)).^2, 3));

%% --- Visualization ---
figure('Name','Painting Differences','Position',[100 100 1200 500]);

subplot(2,2,1);
imagesc(deltaXZ);
axis image off; colorbar; colormap(turbo);
title('ΔXZ');

subplot(2,2,2);
imagesc(deltaXYZ);
axis image off; colorbar; colormap(turbo);
title('ΔXYZ');

subplot(2,2,3);
imagesc(DoGXZ);
axis image off; colorbar; colormap(turbo);
title('DoG XZ');

subplot(2,2,4);
imagesc(DoGXYZ);
axis image off; colorbar; colormap(turbo);
title('DoG XYZ');



%%




