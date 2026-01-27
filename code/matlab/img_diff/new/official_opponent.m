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
        msr_channel = msr_channel - min(msr_channel(:));
        msr_channel = msr_channel / max(msr_channel(:));
        img_msr(:,:,c) = msr_channel;
    end
end



%% painting before vs after
%% --- Setup ---
[h,w,~] = size(painting_before.RGB_lin_img);
min_dim = min(h,w);

% Adaptive scales
sigma_small = 1;
% sigma_large_list = [min_dim/10, min_dim/5, min_dim/3];
sigma_large_list = [min_dim/5, min_dim/3, min_dim/2];
DoG_at = @(X,sigL) imgaussfilt(X,sigma_small) - imgaussfilt(X,sigL);

RGB_before = painting_before.RGB_lin_img;
RGB_after  = painting_after.RGB_lin_img;

%% --- Opponent RGB (De Valois & De Valois, 1988) ---
M_opp = [ 1/sqrt(2)  -1/sqrt(2)   0;
           1/sqrt(6)   1/sqrt(6)  -2/sqrt(6);
           1/sqrt(3)   1/sqrt(3)   1/sqrt(3) ];

Opp_before = reshape(reshape(RGB_before,[],3) * M_opp.', h, w, 3);
Opp_after  = reshape(reshape(RGB_after, [],3) * M_opp.', h, w, 3);

%% --- Δ using chromatic opponent channels (O1, O2) ---
delta_opp = sqrt((Opp_before(:,:,1)-Opp_after(:,:,1)).^2 + ...
                 (Opp_before(:,:,2)-Opp_after(:,:,2)).^2);

%% --- Sum of DoGs across multiple scales ---
DoG_opp = zeros(h,w);
for s = sigma_large_list
    DoG_opp = DoG_opp + abs(DoG_at(Opp_before(:,:,1),s) - DoG_at(Opp_after(:,:,1),s)) + ...
                         abs(DoG_at(Opp_before(:,:,2),s) - DoG_at(Opp_after(:,:,2),s));
end

%% --- Visualization ---
figure('Name','HSI before vs after','Position',[100 100 900 700]);

subplot(1,2,1);
imagesc(delta_opp); axis image off; colorbar; colormap(turbo);
title('Δ - Opponent RGB');

subplot(1,2,2);
imagesc(DoG_opp); axis image off; colorbar; colormap(turbo);
title('sum(DoG) - Opponent RGB');

%% --- Save images with prefix from film file name ---
film_file_prefix = 'film_image'; % replace with actual variable or filename
imwrite(mat2gray(delta_opp), [film_file_prefix '_delta_opp.png']);
imwrite(mat2gray(DoG_opp), [film_file_prefix '_DoG_opp.png']);







%%
%% --- Multi-scale DoG fusion for Opponent RGB chroma (sum only) ---
RGB_film  = film_data.RGB_lin_img;
RGB_hyper = painting_after.RGB_lin_img;

%% --- Prepare 3 MSR cases ---
film_case1  = RGB_film;             
hyper_case1 = RGB_hyper;
film_case2  = MSR(RGB_film, scales); 
hyper_case2 = MSR(RGB_hyper, scales);
film_case3  = MSR(RGB_film, scales); 
hyper_case3 = RGB_hyper;

cases = {
    {'No MSR',        film_case1, hyper_case1};
    {'Both MSR',      film_case2, hyper_case2};
    {'Film only MSR', film_case3, hyper_case3};
};

sigma_small = 1;
[h, w, ~] = size(RGB_film);
min_dim = min(h, w);
% sigma_large_list = [min_dim/10, min_dim/5, min_dim/3];
sigma_large_list = [min_dim/5, min_dim/3, min_dim/2];

nScales = numel(sigma_large_list);

% De Valois orthonormal opponent matrix (O1,O2,L)
M_opp = [ 1/sqrt(2), -1/sqrt(2), 0;
          1/sqrt(6),  1/sqrt(6), -2/sqrt(6);
          1/sqrt(3),  1/sqrt(3),  1/sqrt(3) ];

DoG_at = @(X, sL) imgaussfilt(X, sigma_small) - imgaussfilt(X, sL);

delta_fused_opp = cell(1,3); % sum-fused maps

%% --- Loop over MSR cases ---
for k = 1:3
    rgbF = double(max(cases{k}{2}, 0));
    rgbH = double(max(cases{k}{3}, 0));
    [h,w,~] = size(rgbF);

    % Opponent RGB
    OppF = reshape((reshape(rgbF,[],3) * M_opp.'), h, w, 3);
    OppH = reshape((reshape(rgbH,[],3) * M_opp.'), h, w, 3);

    % Compute DoG maps for each scale
    DoG_maps_opp = zeros(h,w,nScales);
    for si = 1:nScales
        sL = sigma_large_list(si);
        dO1 = DoG_at(OppF(:,:,1), sL) - DoG_at(OppH(:,:,1), sL);
        dO2 = DoG_at(OppF(:,:,2), sL) - DoG_at(OppH(:,:,2), sL);
        DoG_maps_opp(:,:,si) = sqrt(dO1.^2 + dO2.^2);
    end

    % Fuse maps using sum only
    delta_fused_opp{k} = mat2gray(sum(DoG_maps_opp, 3));
end

%% --- Visualization ---
figure('Name','Multi-scale Opponent RGB sum(DoG)','Position',[100 100 1200 400]);
for k = 1:3
    subplot(1,3,k); imagesc(delta_fused_opp{k}); axis image off; colorbar; colormap(turbo);
    title([cases{k}{1}]);
end

%% --- Save images with film filename prefix ---
path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_exp0.mat';
[~, film_prefix, ~] = fileparts(path_film);

for k = 1:3
    fname = [film_prefix '_Opponent_sum_' strrep(cases{k}{1}, ' ', '_') '.png'];
    imwrite(delta_fused_opp{k}, fname);
end








%% Itti
RGB_film  = film_data.RGB_lin_img;
RGB_hyper = painting_after.RGB_lin_img;

%% --- Prepare 3 MSR cases ---
film_case1  = RGB_film;             
hyper_case1 = RGB_hyper;
film_case2  = MSR(RGB_film, scales); 
hyper_case2 = MSR(RGB_hyper, scales);
film_case3  = MSR(RGB_film, scales); 
hyper_case3 = RGB_hyper;

cases = {
    {'No MSR',        film_case1, hyper_case1};
    {'Both MSR',      film_case2, hyper_case2};
    {'Film only MSR', film_case3, hyper_case3};
};

sigma_small = 1;
[h, w, ~] = size(RGB_film);
min_dim = min(h, w);
sigma_large_list = [min_dim/5, min_dim/3, min_dim/2];
nScales = numel(sigma_large_list);

DoG_at = @(X, sL) imgaussfilt(X, sigma_small) - imgaussfilt(X, sL);

delta_fused_opp = cell(1,3); % sum-fused DoG maps
delta_direct_opp = cell(1,3); % Δ maps

%% --- Loop over MSR cases ---
for k = 1:3
    rgbF = double(max(cases{k}{2}, 0));
    rgbH = double(max(cases{k}{3}, 0));
    [h,w,~] = size(rgbF);

    % --- Direct RGB to 2 opponent channels (Itti et al. style) ---
    % O1 = R-G, O2 = B - (R+G)/2
    O1_F = rgbF(:,:,1) - rgbF(:,:,2);
    O2_F = rgbF(:,:,3) - (rgbF(:,:,1)+rgbF(:,:,2))/2;

    O1_H = rgbH(:,:,1) - rgbH(:,:,2);
    O2_H = rgbH(:,:,3) - (rgbH(:,:,1)+rgbH(:,:,2))/2;

    % --- Direct Δ (Euclidean distance in opponent space) ---
    delta_direct_opp{k} = mat2gray(sqrt((O1_F - O1_H).^2 + (O2_F - O2_H).^2));

    % --- Compute DoG maps for each scale ---
    DoG_maps_opp = zeros(h,w,nScales);
    for si = 1:nScales
        sL = sigma_large_list(si);
        dO1 = DoG_at(O1_F, sL) - DoG_at(O1_H, sL);
        dO2 = DoG_at(O2_F, sL) - DoG_at(O2_H, sL);
        DoG_maps_opp(:,:,si) = sqrt(dO1.^2 + dO2.^2);
    end

    % Fuse maps using sum only
    delta_fused_opp{k} = mat2gray(sum(DoG_maps_opp, 3));
end

%% --- Visualization ---
figure('Name','Opponent RGB - Δ and sum(DoG)','Position',[100 100 1200 500]);
for k = 1:3
    subplot(2,3,k); 
    imagesc(delta_direct_opp{k}); axis image off; colorbar; colormap(turbo);
    title([cases{k}{1} ' - Δ']);

    subplot(2,3,3+k); 
    imagesc(delta_fused_opp{k}); axis image off; colorbar; colormap(turbo);
    title([cases{k}{1} ' - sum(DoG)']);
end

%% --- Save images with film filename prefix ---
path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_exp0.mat';
[~, film_prefix, ~] = fileparts(path_film);

for k = 1:3
    imwrite(delta_direct_opp{k}, [film_prefix '_Opponent_delta_' strrep(cases{k}{1}, ' ', '_') '.png']);
    imwrite(delta_fused_opp{k}, [film_prefix '_Opponent_sum_' strrep(cases{k}{1}, ' ', '_') '.png']);
end



%%
%% painting before vs after
%% --- Setup ---
[h,w,~] = size(painting_before.RGB_lin_img);
min_dim = min(h,w);

% Adaptive scales
sigma_small = 1;
sigma_large_list = [min_dim/5, min_dim/3, min_dim/2];
DoG_at = @(X,sigL) imgaussfilt(X,sigma_small) - imgaussfilt(X,sigL);

RGB_before = painting_before.RGB_lin_img;
RGB_after  = painting_after.RGB_lin_img;

%% --- Opponent RGB (direct from RGB, 2 channels) ---
% O1 = R-G, O2 = B-(R+G)/2
O1_before = RGB_before(:,:,1) - RGB_before(:,:,2);
O2_before = RGB_before(:,:,3) - (RGB_before(:,:,1)+RGB_before(:,:,2))/2;

O1_after  = RGB_after(:,:,1) - RGB_after(:,:,2);
O2_after  = RGB_after(:,:,3) - (RGB_after(:,:,1)+RGB_after(:,:,2))/2;

%% --- Δ using chromatic opponent channels (O1, O2) ---
delta_opp = sqrt((O1_before - O1_after).^2 + (O2_before - O2_after).^2);

%% --- Sum of DoGs across multiple scales ---
DoG_opp = zeros(h,w);
for s = sigma_large_list
    DoG_opp = DoG_opp + abs(DoG_at(O1_before,s) - DoG_at(O1_after,s)) + ...
                         abs(DoG_at(O2_before,s) - DoG_at(O2_after,s));
end

%% --- Visualization ---
figure('Name','HSI before vs after','Position',[100 100 900 700]);

subplot(1,2,1);
imagesc(delta_opp); axis image off; colorbar; colormap(turbo);
title('Δ - Opponent RGB');

subplot(1,2,2);
imagesc(DoG_opp); axis image off; colorbar; colormap(turbo);
title('sum(DoG) - Opponent RGB');

%% --- Save images with prefix from film file name ---
path_film = '/Volumes/Study/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_exp0.mat';
[~, film_file_prefix, ~] = fileparts(path_film);

imwrite(mat2gray(delta_opp), [film_file_prefix '_delta_opp.png']);
imwrite(mat2gray(DoG_opp), [film_file_prefix '_DoG_opp.png']);
