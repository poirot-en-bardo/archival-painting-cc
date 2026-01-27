clear; close all;
path_before = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat'; 
path_after = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat'; 
path_before = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat'; 
path_after = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat';



% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat';
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_exp0.mat';
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_underexp.mat';
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/yoda_halogen_fuji_exp0.mat';
% path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/yoda_halogen_fuji_overexp.mat';
path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/yoda_led_kodak_exp0.mat';

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


%% Step 1: Multi-Scale Retinex (MSR)
function img_msr = MSR(img, scales)
    img_msr = zeros(size(img));
    img = max(img,1e-6); % avoid log(0)
    weight = 1/length(scales);
    for c = 1:3
        channel = img(:,:,c);
        msr_channel = zeros(size(channel));
        for s = 1:length(scales)
            sigma = scales(s);
            blur = imgaussfilt(channel, sigma);
            SSR = log(channel) - log(blur);
            msr_channel = msr_channel + weight*SSR;
        end
        % msr_channel = msr_channel - min(msr_channel(:));
        % msr_channel = msr_channel / max(msr_channel(:));
        eps_val = 1e-8;
        msr_channel = (msr_channel - min(msr_channel(:))) / (max(msr_channel(:)) - min(msr_channel(:)) + eps_val);

        img_msr(:,:,c) = msr_channel;
    end
end


%% --- Setup ---
RGB_film   = film_data.RGB_lin_img;
RGB_hyper  = painting_after.RGB_lin_img;
RGB_before = painting_before.RGB_lin_img;
RGB_after  = painting_after.RGB_lin_img;

[h,w,~] = size(RGB_film);
min_dim = min(h,w);
sigma_small = 1;
% sigma_large_list = [min_dim/5, min_dim/3, min_dim/2]; 
sigma_large_list = [min_dim/10, min_dim/5, min_dim/3]; 
nScales = numel(sigma_large_list);

DoG_at = @(X, sL) imgaussfilt(X, sigma_small) - imgaussfilt(X, sL);

%% --- Prepare 3 MSR cases (film vs HSI) ---
scales = [15,80,250];
film_case1  = RGB_film;
hyper_case1 = RGB_hyper;
film_case2  = MSR(RGB_film, scales);
hyper_case2 = MSR(RGB_hyper, scales);
film_case3  = MSR(RGB_film, scales);
hyper_case3 = RGB_hyper;

cases_film_hsi = {
    {'No MSR',        film_case1, hyper_case1};
    {'Both MSR',      film_case2, hyper_case2};
    {'Film only MSR', film_case3, hyper_case3};
};

%% --- Initialize storage ---
delta_YCbCr_film = cell(1,3);
DoG_fused_YCbCr_film = cell(1,3);

%% --- Compute Δ and mean-of-DoG for film vs HSI ---
for k = 1:3
    rgbF = cases_film_hsi{k}{2};
    rgbH = cases_film_hsi{k}{3};
    
    % YCbCr
    YCbCr_F = rgb2ycbcr(rgbF);
    YCbCr_H = rgb2ycbcr(rgbH);

    % Δ chromatic map
    delta_YCbCr_film{k} = mat2gray(sqrt( ...
        (YCbCr_F(:,:,2)-YCbCr_H(:,:,2)).^2 + ...
        (YCbCr_F(:,:,3)-YCbCr_H(:,:,3)).^2 ));

    % mean-of-DoG map
    DoG_maps = zeros(h,w,nScales);
    for si = 1:nScales
        sL = sigma_large_list(si);
        dC1 = DoG_at(YCbCr_F(:,:,2), sL) - DoG_at(YCbCr_H(:,:,2), sL);
        dC2 = DoG_at(YCbCr_F(:,:,3), sL) - DoG_at(YCbCr_H(:,:,3), sL);
        DoG_maps(:,:,si) = sqrt(dC1.^2 + dC2.^2);
    end
    DoG_fused_YCbCr_film{k} = mat2gray(mean(DoG_maps,3));
end

%% --- Compute Δ and mean-of-DoG for before vs after ---
YCbCr_B = rgb2ycbcr(RGB_before);
YCbCr_A = rgb2ycbcr(RGB_after);

% Δ chromatic map
delta_YCbCr_before_after = mat2gray(sqrt( ...
    (YCbCr_B(:,:,2)-YCbCr_A(:,:,2)).^2 + ...
    (YCbCr_B(:,:,3)-YCbCr_A(:,:,3)).^2 ));

% mean-of-DoG map
DoG_maps = zeros(h,w,nScales);
for si = 1:nScales
    sL = sigma_large_list(si);
    dC1 = DoG_at(YCbCr_B(:,:,2), sL) - DoG_at(YCbCr_A(:,:,2), sL);
    dC2 = DoG_at(YCbCr_B(:,:,3), sL) - DoG_at(YCbCr_A(:,:,3), sL);
    DoG_maps(:,:,si) = sqrt(dC1.^2 + dC2.^2);
end
DoG_YCbCr_before_after = mat2gray(mean(DoG_maps,3));

%% --- Visualization ---

% Film vs HSI
figure('Name','Film vs HSI','Position',[50 50 1400 700]);
for k = 1:3
    subplot(2,3,k);
    imagesc(delta_YCbCr_film{k}); axis image off; colorbar; colormap(turbo);
    title([cases_film_hsi{k}{1} ' Δ']);
    
    subplot(2,3,3+k);
    imagesc(DoG_fused_YCbCr_film{k}); axis image off; colorbar; colormap(turbo);
    title([cases_film_hsi{k}{1} ' mean-of-DoG']);
end

% Before vs After
figure('Name','Before vs After','Position',[50 50 1200 500]);
subplot(1,2,1);
imagesc(delta_YCbCr_before_after); axis image off; colorbar; colormap(turbo);
title('Δ Before vs After');

subplot(1,2,2);
imagesc(DoG_YCbCr_before_after); axis image off; colorbar; colormap(turbo);
title('mean-of-DoG Before vs After');





%% --- Binary masks based on percentile ---

percentile_val = 85; % 85th percentile

% --- Before vs After ---
thresh_BA = prctile(DoG_YCbCr_before_after(:), percentile_val);
mask_BA = DoG_YCbCr_before_after >= thresh_BA;

figure('Name','Binary Mask: Before vs After','Position',[100 100 600 500]);
imshow(mask_BA);
title(['Binary mask Before vs After (DoG >= ' num2str(percentile_val) 'th percentile)']);

% --- Film vs After (No MSR only) ---
DoG_film_noMSR = DoG_fused_YCbCr_film{1};
thresh_FA = prctile(DoG_film_noMSR(:), percentile_val);
mask_FA = DoG_film_noMSR >= thresh_FA;

figure('Name','Binary Mask: Film vs After (No MSR)','Position',[100 100 600 500]);
imshow(mask_FA);
title(['Binary mask Film vs After (No MSR, DoG >= ' num2str(percentile_val) 'th percentile)']);

