%% Inputs
% film_rgb: film photograph (double, [0,1])
% hyper_rgb: hyperspectral-derived RGB image (double, [0,1])

path_before = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat';
path_after  = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat';
path_before = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat';
path_after  = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat';

path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat';
% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_exp0.mat';
% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_underexp.mat';
path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/yoda_halogen_fuji_exp0.mat';
% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/yoda_halogen_fuji_overexp.mat';
% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/yoda_led_kodak_exp0.mat';




painting_before = load(path_before);
painting_after  = load(path_after);
film_data       = load(path_film);

img_before = painting_before.RGB_img;
img_after   = painting_after.RGB_img;
film_rgb    = film_data.RGB_img;
hyper_rgb = img_after;

% ground truth
% Lab_before = painting_before.Lab_img;
% Lab_after  = painting_after.Lab_img;
% dE_map_direct = reshape(deltaE2000(reshape(Lab_before,[],3), reshape(Lab_after,[],3)), size(Lab_before,1), size(Lab_before,2));


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
        % Normalize to 0-100 for Lab scale consistency
        msr_channel = msr_channel - min(msr_channel(:));
        msr_channel = msr_channel / max(msr_channel(:));
        img_msr(:,:,c) = msr_channel;
    end
end

scales = [15, 80, 250];
film_msr  = MSR(film_rgb, scales);
hyper_msr = MSR(hyper_rgb, scales);

%% ======================== Step 2: Convert to Lab ========================
isGamma = true;
XYZ_film  = prophoto2xyz(film_msr, isGamma);
XYZ_hyper = prophoto2xyz(hyper_msr, isGamma);

film_lab  = xyz2lab_custom(reshape(XYZ_film,[],3));
hyper_lab = xyz2lab_custom(reshape(XYZ_hyper,[],3));

% Reshape back to image
film_lab  = reshape(film_lab, size(film_rgb));
hyper_lab = reshape(hyper_lab, size(hyper_rgb));

L_film   = film_lab(:,:,1); 
L_hyper  = hyper_lab(:,:,1);


%% ======================== Step 3: Luminance Normalization ========================
p1 = prctile(L_film(:),1); p99 = prctile(L_film(:),99);
L_film_norm  = 100*(L_film - p1)/(p99-p1); 
L_film_norm  = min(max(L_film_norm,0),100);

p1 = prctile(L_hyper(:),1); p99 = prctile(L_hyper(:),99);
L_hyper_norm = 100*(L_hyper - p1)/(p99-p1); 
L_hyper_norm = min(max(L_hyper_norm,0),100);

%% ======================== Step 4: DoG Filtering ========================
sigma_small = 3; sigma_large = 15; % larger sigma for Lab L range
DoG_film  = imgaussfilt(L_film_norm, sigma_small) - imgaussfilt(L_film_norm, sigma_large);
DoG_hyper = imgaussfilt(L_hyper_norm, sigma_small) - imgaussfilt(L_hyper_norm, sigma_large);
DoG_diff  = abs(DoG_film - DoG_hyper);

%% ======================== Step 5: Gradient-based Difference ========================
[Gx_film, Gy_film]   = imgradientxy(L_film_norm);
[Gx_hyper, Gy_hyper] = imgradientxy(L_hyper_norm);
grad_diff = sqrt((Gx_film - Gx_hyper).^2 + (Gy_film - Gy_hyper).^2);

% Combine DoG and gradient differences
alpha = 0.5;
change_map = alpha*DoG_diff + (1-alpha)*grad_diff;

%% ======================== Step 6: Threshold + Cleanup ========================
threshold = 2; % adjust based on Lab/DoG scale (0-100)
change_mask = change_map > threshold;

change_mask = bwareaopen(change_mask,50);          % remove small noise
change_mask = imclose(change_mask, strel('disk',3)); % smooth edges

%% ======================== Step 7: Visualization ========================
figure; imshow(film_rgb); title('Detected Material Changes Overlay');
hold on;
h = imshow(cat(3, ones(size(change_mask)), zeros(size(change_mask)), zeros(size(change_mask))));
set(h,'AlphaData',0.5*change_mask);

%% ======================== Step 8: Optional: Intermediate Maps ========================
figure;
subplot(2,2,1); imshow(L_film_norm,[]); title('Normalized Film L');
subplot(2,2,2); imshow(L_hyper_norm,[]); title('Normalized Hyper L');
subplot(2,2,3); imshow(DoG_diff,[]); title('DoG Difference');
subplot(2,2,4); imshow(grad_diff,[]); title('Gradient Difference');

%%
% Example: display DoG difference, clamped to 0–200
figure;
imagesc(DoG_diff, [0 5]);   % display values from 0 to 200
axis image; axis off; colorbar;
title('DoG Difference (0-200 scale)');

% Gradient difference, clamped similarly
figure;
imagesc(grad_diff, [0 10]);
axis image; axis off; colorbar;
title('Gradient Difference (0-200 scale)');
