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

%%
% ground truth
Lab_before = painting_before.Lab_img;
Lab_after  = painting_after.Lab_img;
Lab_film = film_data.Lab_img;
% dE_map_direct = reshape(deltaE2000(reshape(Lab_before,[],3), reshape(Lab_after,[],3)), size(Lab_before,1), size(Lab_before,2));


%% Inputs
% film_rgb: film photograph (double, [0,1])
% hyper_rgb: hyperspectral-derived RGB image (double, [0,1])

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

scales = [15, 80, 250];
film_msr = MSR(film_rgb, scales);
hyper_msr = MSR(hyper_rgb, scales);

% Intermediate visualization
figure;
subplot(2,2,1); imshow(film_rgb); title('Original Film');
subplot(2,2,2); imshow(film_msr); title('Film after MSR');
subplot(2,2,3); imshow(hyper_rgb); title('Original Hyper');
subplot(2,2,4); imshow(hyper_msr); title('Hyper after MSR');

% Intermediate DoG difference after MSR
DoG_temp = abs(imgaussfilt(rgb2gray(film_msr),1)-imgaussfilt(rgb2gray(hyper_msr),1));
figure; imshow(DoG_temp,[]); title('Intermediate DoG Diff (after MSR)');

%% Step 2: Convert to Lab
isGamma = true;
XYZ_film_msr  = prophoto2xyz(film_msr, isGamma);
XYZ_hyper_msr = prophoto2xyz(hyper_msr, isGamma);

film_lab_msr  = xyz2lab_custom(reshape(XYZ_film_msr,[],3));
hyper_lab_msr = xyz2lab_custom(reshape(XYZ_hyper_msr,[],3));

% Reshape back to image
film_lab_msr  = reshape(film_lab_msr, size(film_rgb));
hyper_lab_msr = reshape(hyper_lab_msr, size(hyper_rgb));

% Extract channels
L_film_msr   = film_lab_msr(:,:,1); 
a_film_msr   = film_lab_msr(:,:,2); 
b_film_msr   = film_lab_msr(:,:,3);

L_hyper_msr  = hyper_lab_msr(:,:,1); 
a_hyper_msr  = hyper_lab_msr(:,:,2); 
b_hyper_msr  = hyper_lab_msr(:,:,3);

% Visualization
figure;
subplot(2,3,1); imshow(L_film_msr,[]); title('Film L (after MSR)');
subplot(2,3,2); imshow(a_film_msr,[]); title('Film a (after MSR)');
subplot(2,3,3); imshow(b_film_msr,[]); title('Film b (after MSR)');
subplot(2,3,4); imshow(L_hyper_msr,[]); title('Hyper L (after MSR)');
subplot(2,3,5); imshow(a_hyper_msr,[]); title('Hyper a (after MSR)');
subplot(2,3,6); imshow(b_hyper_msr,[]); title('Hyper b (after MSR)');

%%
% Extract Lab channels before MSR
L_film   = Lab_film(:,:,1); 
a_film   = Lab_film(:,:,2); 
b_film   = Lab_film(:,:,3);

L_hyper  = Lab_after(:,:,1); 
a_hyper  = Lab_after(:,:,2); 
b_hyper  = Lab_after(:,:,3);

% Visualization of Lab channels before MSR
figure;
subplot(2,3,1); imshow(L_film, []); title('Film L (before MSR)');
subplot(2,3,2); imshow(a_film, []); title('Film a (before MSR)');
subplot(2,3,3); imshow(b_film, []); title('Film b (before MSR)');
subplot(2,3,4); imshow(L_hyper, []); title('Hyper L (before MSR)');
subplot(2,3,5); imshow(a_hyper, []); title('Hyper a (before MSR)');
subplot(2,3,6); imshow(b_hyper, []); title('Hyper b (before MSR)');

%%
% Intermediate change map in L channel
L_diff1 = abs(L_film - L_hyper);
L_diff1 = abs(a_film - a_hyper);
% L_diff1 = abs(b_film - b_hyper);

figure; imshow(L_diff1,[]); title('Intermediate L Difference (Lab)');

%% Step 3: Percentile-based normalization (exposure compensation)
p1 = prctile(L_film(:),1); p99 = prctile(L_film(:),99);
L_film_norm = (L_film - p1)/(p99-p1); L_film_norm = min(max(L_film_norm,0),1);

p1 = prctile(L_hyper(:),1); p99 = prctile(L_hyper(:),99);
L_hyper_norm = (L_hyper - p1)/(p99-p1); L_hyper_norm = min(max(L_hyper_norm,0),1);

% Visualize normalized L and intermediate difference
figure;
subplot(1,3,1); imshow(L_film_norm); title('Normalized Film L');
subplot(1,3,2); imshow(L_hyper_norm); title('Normalized Hyper L');
subplot(1,3,3); imshow(abs(L_film_norm-L_hyper_norm),[]); title('L Difference (Normalized)');

%% Step 4: DoG filtering
sigma_small = 1; sigma_large = 5;
DoG_film = imgaussfilt(L_film_norm,sigma_small) - imgaussfilt(L_film_norm,sigma_large);
DoG_hyper = imgaussfilt(L_hyper_norm,sigma_small) - imgaussfilt(L_hyper_norm,sigma_large);

DoG_diff = abs(DoG_film - DoG_hyper);

figure;
subplot(1,3,1); imshow(DoG_film,[]); title('DoG Film');
subplot(1,3,2); imshow(DoG_hyper,[]); title('DoG Hyper');
subplot(1,3,3); imshow(DoG_diff,[]); title('DoG Difference');

%% Step 5: Gradient-based difference
[Gx_film, Gy_film] = imgradientxy(L_film_norm);
[Gx_hyper, Gy_hyper] = imgradientxy(L_hyper_norm);
grad_diff = sqrt((Gx_film - Gx_hyper).^2 + (Gy_film - Gy_hyper).^2);

figure;
subplot(1,3,1); imshow(Gx_film,[]); title('Gradient X Film');
subplot(1,3,2); imshow(Gy_film,[]); title('Gradient Y Film');
subplot(1,3,3); imshow(grad_diff,[]); title('Gradient Difference');

%% Step 6: Combine DoG and gradient differences
alpha = 0.5;
change_map = alpha*DoG_diff + (1-alpha)*grad_diff;

figure; imshow(change_map,[]); title('Combined Change Map');

%% Step 7: Threshold and cleanup
threshold = 0.05;
change_mask = change_map > threshold;
change_mask = bwareaopen(change_mask,50);
change_mask = imclose(change_mask, strel('disk',3));

figure; imshow(change_mask); title('Final Binary Change Mask');

%% Step 8: Overlay on original film image
figure; imshow(film_rgb); title('Detected Material Changes Overlay');
hold on;
h = imshow(cat(3, ones(size(change_mask)), zeros(size(change_mask)), zeros(size(change_mask))));
set(h,'AlphaData',0.5*change_mask);

%%
% Before MSR
ab_film_nomsr  = cat(3, a_film, b_film);
ab_hyper_nomsr = cat(3, a_hyper, b_hyper);

% After MSR
ab_film_msr  = cat(3, a_film_msr, b_film_msr);
ab_hyper_msr = cat(3, a_hyper_msr, b_hyper_msr);

% Before MSR
delta_ab_nomsr = sqrt( (ab_film_nomsr(:,:,1) - ab_hyper_nomsr(:,:,1)).^2 + ...
                       (ab_film_nomsr(:,:,2) - ab_hyper_nomsr(:,:,2)).^2 );

% After MSR
delta_ab_msr   = sqrt( (ab_film_msr(:,:,1) - ab_hyper_msr(:,:,1)).^2 + ...
                       (ab_film_msr(:,:,2) - ab_hyper_msr(:,:,2)).^2 );
%%
threshold = 15; % adjust based on perceptual sensitivity

change_mask_ab_nomsr = delta_ab_nomsr > threshold;
change_mask_ab_msr   = delta_ab_msr > threshold;

% Visualize
figure;
subplot(1,2,1); imshow(change_mask_ab_nomsr); title('Material Change Mask (a*b*, before MSR)');
subplot(1,2,2); imshow(change_mask_ab_msr); title('Material Change Mask (a*b*, after MSR)');
