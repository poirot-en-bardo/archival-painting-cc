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
        msr_channel = msr_channel - min(msr_channel(:));
        msr_channel = msr_channel / max(msr_channel(:));
        img_msr(:,:,c) = msr_channel;
    end
end





RGB_film  = film_data.RGB_lin_img;
RGB_hyper = painting_after.RGB_lin_img;

%% --- Prepare 3 MSR cases ---
scales = [15, 80, 250];  % MSR scales
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
nScales = numel(sigma_lar]= ...
    ge_list);

DoG_at = @(X, sL) imgaussfilt(X, sigma_small) - imgaussfilt(X, sL);

delta_fused_opp = cell(1,3); % sum-fused DoG maps
delta_direct_opp = cell(1,3); % Δ maps

delta_fused_lab = cell(1,3);      % sum-fused DoG maps
delta_direct_lab = cell(1,3);     % Δ maps

for k = 1:3
    rgbF = double(max(cases{k}{2}, 0));
    rgbH = double(max(cases{k}{3}, 0));
    [h,w,~] = size(rgbF);

    % --- Convert ProPhoto RGB to XYZ ---
    XYZF = prophoto2xyz(rgbF, false);
    XYZH = prophoto2xyz(rgbH, false);

    % --- Reshape to Nx3 for Lab conversion ---
    XYZF_reshaped = reshape(XYZF, [], 3);
    XYZH_reshaped = reshape(XYZH, [], 3);

    labF_reshaped = xyz2lab_custom(XYZF_reshaped);  % Nx3
    labH_reshaped = xyz2lab_custom(XYZH_reshaped);

    % --- Reshape back to HxWx3 ---
    labF = reshape(labF_reshaped, h, w, 3);
    labH = reshape(labH_reshaped, h, w, 3);

    % --- a & b channels ---
    a_F = labF(:,:,2);
    b_F = labF(:,:,3);

    a_H = labH(:,:,2);
    b_H = labH(:,:,3);

    % --- Direct Δ (Euclidean distance in a-b space) ---
    delta_direct_lab{k} = mat2gray(sqrt((a_F - a_H).^2 + (b_F - b_H).^2));

    % --- Compute DoG maps for each scale ---
    DoG_maps_lab = zeros(h,w,nScales);
    for si = 1:nScales
        sL = sigma_large_list(si);
        da = DoG_at(a_F, sL) - DoG_at(a_H, sL);
        db = DoG_at(b_F, sL) - DoG_at(b_H, sL);
        DoG_maps_lab(:,:,si) = sqrt(da.^2 + db.^2);
    end

    % --- Sum fusion only ---
    delta_fused_lab{k} = mat2gray(sum(DoG_maps_lab, 3));
end

%% --- Visualization ---
figure('Name','Lab a-b - Δ and sum(DoG)','Position',[100 100 1200 500]);
for k = 1:3
    subplot(2,3,k); 
    imagesc(delta_direct_lab{k}); axis image off; colorbar; colormap(parula);
    title([cases{k}{1} ' - Δ']);

    subplot(2,3,3+k); 
    imagesc(delta_fused_lab{k}); axis image off; colorbar; colormap(parula);
    title([cases{k}{1} ' - sum(DoG)']);
end
