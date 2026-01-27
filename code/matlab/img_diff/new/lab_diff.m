path_before = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat';
path_after  = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat';
path_before = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat';
path_after  = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat';

% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat';
% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_exp0.mat';
% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_underexp.mat';
% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/yoda_halogen_fuji_exp0.mat';
% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/yoda_halogen_fuji_overexp.mat';
path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/yoda_led_kodak_exp0.mat';

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

scales = [1, 20, 200];
film_msr = MSR(film_rgb, scales);
hyper_msr = MSR(hyper_rgb, scales);

% Intermediate visualization
figure;
subplot(2,2,1); imshow(film_rgb); title('Original Film');
subplot(2,2,2); imshow(film_msr); title('Film after MSR');
subplot(2,2,3); imshow(hyper_rgb); title('Original Hyper');
subplot(2,2,4); imshow(hyper_msr); title('Hyper after MSR');



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
subplot(2,3,6); imshow(b_hyper, []); title('Hy per b (before MSR)');

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
threshold = 5; % adjust based on perceptual sensitivity

change_mask_ab_nomsr = delta_ab_nomsr > threshold;
change_mask_ab_msr   = delta_ab_msr > threshold;

%%
p = 85; % percentile
threshold_nomsr = prctile(delta_ab_nomsr(:), p);
threshold_msr   = prctile(delta_ab_msr(:), p);

change_mask_ab_nomsr = delta_ab_nomsr > threshold_nomsr;
change_mask_ab_msr   = delta_ab_msr   > threshold_msr;
% Visualize
% figure;
% subplot(1,2,1); imshow(change_mask_ab_nomsr); title('Material Change Mask (a*b*, before MSR)');
% subplot(1,2,2); imshow(change_mask_ab_msr); title('Material Change Mask (a*b*, after MSR)');

%%
k = 1; % multiplier
threshold_nomsr = mean(delta_ab_nomsr(:)) + k*std(delta_ab_nomsr(:));
threshold_msr   = mean(delta_ab_msr(:))   + k*std(delta_ab_msr(:));

change_mask_ab_nomsr = delta_ab_nomsr > threshold_nomsr;
change_mask_ab_msr   = delta_ab_msr   > threshold_msr;



% Visualize
% figure;
% subplot(1,2,1); imshow(change_mask_ab_nomsr); title('Material Change Mask (a*b*, before MSR)');
% subplot(1,2,2); imshow(change_mask_ab_msr); title('Material Change Mask (a*b*, after MSR)');

%% dog on Lab
% Extract a* and b* channels

% Before MSR
a_film_nomsr  = Lab_film(:,:,2);
b_film_nomsr  = Lab_film(:,:,3);

a_hyper_nomsr = Lab_after(:,:,2);
b_hyper_nomsr = Lab_after(:,:,3);

% After MSR
a_film_msr  = film_lab_msr(:,:,2);
b_film_msr  = film_lab_msr(:,:,3);

a_hyper_msr = hyper_lab_msr(:,:,2);
b_hyper_msr = hyper_lab_msr(:,:,3);

% Apply DoG filtering
sigma_small = 1;
sigma_large = 500;

% Before MSR
DoG_a_film_nomsr  = imgaussfilt(a_film_nomsr, sigma_small) - imgaussfilt(a_film_nomsr, sigma_large);
DoG_b_film_nomsr  = imgaussfilt(b_film_nomsr, sigma_small) - imgaussfilt(b_film_nomsr, sigma_large);
DoG_a_hyper_nomsr = imgaussfilt(a_hyper_nomsr, sigma_small) - imgaussfilt(a_hyper_nomsr, sigma_large);
DoG_b_hyper_nomsr = imgaussfilt(b_hyper_nomsr, sigma_small) - imgaussfilt(b_hyper_nomsr, sigma_large);

% After MSR
DoG_a_film_msr  = imgaussfilt(a_film_msr, sigma_small) - imgaussfilt(a_film_msr, sigma_large);
DoG_b_film_msr  = imgaussfilt(b_film_msr, sigma_small) - imgaussfilt(b_film_msr, sigma_large);
DoG_a_hyper_msr = imgaussfilt(a_hyper_msr, sigma_small) - imgaussfilt(a_hyper_msr, sigma_large);
DoG_b_hyper_msr = imgaussfilt(b_hyper_msr, sigma_small) - imgaussfilt(b_hyper_msr, sigma_large);

% Compute absolute differences and combine channels
DoG_diff_ab_nomsr = sqrt( (DoG_a_film_nomsr - DoG_a_hyper_nomsr).^2 + ...
                           (DoG_b_film_nomsr - DoG_b_hyper_nomsr).^2 );

DoG_diff_ab_msr   = sqrt( (DoG_a_film_msr - DoG_a_hyper_msr).^2 + ...
                           (DoG_b_film_msr - DoG_b_hyper_msr).^2 );

% Normalize for visualization
DoG_diff_ab_nomsr = DoG_diff_ab_nomsr - min(DoG_diff_ab_nomsr(:));
DoG_diff_ab_nomsr = DoG_diff_ab_nomsr / max(DoG_diff_ab_nomsr(:));

DoG_diff_ab_msr = DoG_diff_ab_msr - min(DoG_diff_ab_msr(:));
DoG_diff_ab_msr = DoG_diff_ab_msr / max(DoG_diff_ab_msr(:));

%%
% Display side-by-side
figure;
subplot(1,2,1); imagesc(DoG_diff_ab_nomsr, [0 0.5]); axis image; axis off;
colormap(turbo); colorbar;
title('DoG on a*b* (before MSR)');

subplot(1,2,2); imagesc(DoG_diff_ab_msr); axis image; axis off;
colormap(turbo); colorbar;
title('DoG on a*b* (after MSR)');

%%


%% Assume we already have:
% film_lab_msr  -> [M x N x 3]
% hyper_lab_msr -> [M x N x 3]

[M,N,~] = size(film_lab_msr);

% Reshape to columns (Nx3) for deltaE2000
LAB_film_cols  = reshape(film_lab_msr, [], 3);  % (M*N) x 3
LAB_hyper_cols = reshape(hyper_lab_msr, [], 3);

% Optional: use weighting factors if needed, e.g., K = [1 1 1];
K = [1 1 1];

% Compute Delta E
DE_cols = deltaE2000(LAB_hyper_cols, LAB_film_cols, K);  % reference = hyperspectral

% Reshape back to image
DE_img = reshape(DE_cols, M, N);

% Visualization
figure;
imagesc(DE_img, [0 20]); axis image off; colorbar;
colormap(jet);
title('\DeltaE_{00} (Film vs Hyperspectral, after MSR)');

%%
% Δab before MSR
delta_ab_nomsr = sqrt((a_film - a_hyper).^2 + (b_film - b_hyper).^2);

% Δab after MSR
delta_ab_msr   = sqrt((a_film_msr - a_hyper_msr).^2 + (b_film_msr - b_hyper_msr).^2);

% Visualization
figure;
subplot(1,2,1); imagesc(delta_ab_nomsr, [0 60]); axis image off; colorbar;
colormap(jet); title('\Delta a*b* (before MSR)');

subplot(1,2,2); imagesc(delta_ab_msr); axis image off; colorbar;
colormap(jet); title('\Delta a*b* (after MSR)');


%%

%% Compute delta_CbCr (before and after MSR)

% --- Before MSR ---
film_ycbcr_nomsr  = rgb2ycbcr(film_rgb);
hyper_ycbcr_nomsr = rgb2ycbcr(hyper_rgb);

Cb_film_nomsr  = double(film_ycbcr_nomsr(:,:,2));
Cr_film_nomsr  = double(film_ycbcr_nomsr(:,:,3));

Cb_hyper_nomsr = double(hyper_ycbcr_nomsr(:,:,2));
Cr_hyper_nomsr = double(hyper_ycbcr_nomsr(:,:,3));

delta_CbCr_nomsr = sqrt((Cb_film_nomsr - Cb_hyper_nomsr).^2 + ...
                        (Cr_film_nomsr - Cr_hyper_nomsr).^2);

% --- After MSR ---
film_ycbcr_msr  = rgb2ycbcr(film_msr);
hyper_ycbcr_msr = rgb2ycbcr(hyper_rgb); %change

Cb_film_msr  = double(film_ycbcr_msr(:,:,2));
Cr_film_msr  = double(film_ycbcr_msr(:,:,3));

Cb_hyper_msr = double(hyper_ycbcr_msr(:,:,2));
Cr_hyper_msr = double(hyper_ycbcr_msr(:,:,3));

delta_CbCr_msr = sqrt((Cb_film_msr - Cb_hyper_msr).^2 + ...
                      (Cr_film_msr - Cr_hyper_msr).^2);

% Visualization
figure;
subplot(1,2,1);
imagesc(delta_CbCr_nomsr); axis image off; colorbar;
colormap(jet); title('\Delta CbCr (before MSR)');

subplot(1,2,2);
imagesc(delta_CbCr_msr); axis image off; colorbar;
colormap(jet); title('\Delta CbCr (after MSR)');

%%
%% Visualize YCbCr channels (before and after MSR)
% --- Before MSR ---
film_ycbcr_nomsr  = rgb2ycbcr(film_rgb);
hyper_ycbcr_nomsr = rgb2ycbcr(hyper_rgb);

Y_film_nomsr  = film_ycbcr_nomsr(:,:,1);
Cb_film_nomsr = film_ycbcr_nomsr(:,:,2);
Cr_film_nomsr = film_ycbcr_nomsr(:,:,3);

Y_hyper_nomsr  = hyper_ycbcr_nomsr(:,:,1);
Cb_hyper_nomsr = hyper_ycbcr_nomsr(:,:,2);
Cr_hyper_nomsr = hyper_ycbcr_nomsr(:,:,3);

figure;
subplot(2,3,1); imshow(Y_film_nomsr,[]); title('Film Y (before MSR)');
subplot(2,3,2); imshow(Cb_film_nomsr,[]); title('Film Cb (before MSR)');
subplot(2,3,3); imshow(Cr_film_nomsr,[]); title('Film Cr (before MSR)');
subplot(2,3,4); imshow(Y_hyper_nomsr,[]); title('Hyper Y (before MSR)');
subplot(2,3,5); imshow(Cb_hyper_nomsr,[]); title('Hyper Cb (before MSR)');
subplot(2,3,6); imshow(Cr_hyper_nomsr,[]); title('Hyper Cr (before MSR)');

% --- After MSR ---
film_ycbcr_msr  = rgb2ycbcr(film_msr);
hyper_ycbcr_msr = rgb2ycbcr(hyper_msr);

Y_film_msr  = film_ycbcr_msr(:,:,1);
Cb_film_msr = film_ycbcr_msr(:,:,2);
Cr_film_msr = film_ycbcr_msr(:,:,3);

Y_hyper_msr  = hyper_ycbcr_msr(:,:,1);
Cb_hyper_msr = hyper_ycbcr_msr(:,:,2);
Cr_hyper_msr = hyper_ycbcr_msr(:,:,3);

figure;
subplot(2,3,1); imshow(Y_film_msr,[]); title('Film Y (after MSR)');
subplot(2,3,2); imshow(Cb_film_msr,[]); title('Film Cb (after MSR)');
subplot(2,3,3); imshow(Cr_film_msr,[]); title('Film Cr (after MSR)');
subplot(2,3,4); imshow(Y_hyper_msr,[]); title('Hyper Y (after MSR)');
subplot(2,3,5); imshow(Cb_hyper_msr,[]); title('Hyper Cb (after MSR)');
subplot(2,3,6); imshow(Cr_hyper_msr,[]); title('Hyper Cr (after MSR)');


%%
%% DoG on CbCr channels (before and after MSR)

sigma_small = 1;
sigma_large = 300;

% --- Before MSR ---
DoG_Cb_film_nomsr  = imgaussfilt(Cb_film_nomsr, sigma_small) - imgaussfilt(Cb_film_nomsr, sigma_large);
DoG_Cr_film_nomsr  = imgaussfilt(Cr_film_nomsr, sigma_small) - imgaussfilt(Cr_film_nomsr, sigma_large);
DoG_Cb_hyper_nomsr = imgaussfilt(Cb_hyper_nomsr, sigma_small) - imgaussfilt(Cb_hyper_nomsr, sigma_large);
DoG_Cr_hyper_nomsr = imgaussfilt(Cr_hyper_nomsr, sigma_small) - imgaussfilt(Cr_hyper_nomsr, sigma_large);

DoG_diff_CbCr_nomsr = sqrt( (DoG_Cb_film_nomsr - DoG_Cb_hyper_nomsr).^2 + ...
                            (DoG_Cr_film_nomsr - DoG_Cr_hyper_nomsr).^2 );

% --- After MSR ---
DoG_Cb_film_msr  = imgaussfilt(Cb_film_msr, sigma_small) - imgaussfilt(Cb_film_msr, sigma_large);
DoG_Cr_film_msr  = imgaussfilt(Cr_film_msr, sigma_small) - imgaussfilt(Cr_film_msr, sigma_large);
DoG_Cb_hyper_msr = imgaussfilt(Cb_hyper_msr, sigma_small) - imgaussfilt(Cb_hyper_msr, sigma_large);
DoG_Cr_hyper_msr = imgaussfilt(Cr_hyper_msr, sigma_small) - imgaussfilt(Cr_hyper_msr, sigma_large);

DoG_diff_CbCr_msr = sqrt( (DoG_Cb_film_msr - DoG_Cb_hyper_msr).^2 + ...
                          (DoG_Cr_film_msr - DoG_Cr_hyper_msr).^2 );

% Normalize for visualization
DoG_diff_CbCr_nomsr = mat2gray(DoG_diff_CbCr_nomsr);
DoG_diff_CbCr_msr   = mat2gray(DoG_diff_CbCr_msr);

% --- Visualization ---
figure;
subplot(1,2,1);
imagesc(DoG_diff_CbCr_nomsr, [0 1]); axis image off; colormap(turbo); colorbar;
title('DoG on CbCr (before MSR)');

subplot(1,2,2);
imagesc(DoG_diff_CbCr_msr, [0 1]); axis image off; colormap(turbo); colorbar;
title('DoG on CbCr (after MSR)');

%%







%% 🔁 Cross comparison: Film after MSR vs Hyper before MSR

% Define your "reference" pair
film_ref   = film_msr;     % corrected film
hyper_ref  = hyper_rgb;    % uncorrected hyperspectral RGB

% Convert both to Lab
isGamma = true;
XYZ_film_ref  = prophoto2xyz(film_ref, isGamma);
XYZ_hyper_ref = prophoto2xyz(hyper_ref, isGamma);

film_lab_ref  = xyz2lab_custom(reshape(XYZ_film_ref,[],3));
hyper_lab_ref = xyz2lab_custom(reshape(XYZ_hyper_ref,[],3));

film_lab_ref  = reshape(film_lab_ref, size(film_ref));
hyper_lab_ref = reshape(hyper_lab_ref, size(hyper_ref));

% Extract channels
a_film_ref = film_lab_ref(:,:,2);
b_film_ref = film_lab_ref(:,:,3);
a_hyper_ref = hyper_lab_ref(:,:,2);
b_hyper_ref = hyper_lab_ref(:,:,3);

% 3️⃣ ΔCbCr
film_ycbcr_ref  = rgb2ycbcr(film_ref);
hyper_ycbcr_ref = rgb2ycbcr(hyper_ref);

Cb_film_ref = double(film_ycbcr_ref(:,:,2));
Cr_film_ref = double(film_ycbcr_ref(:,:,3));
Cb_hyper_ref = double(hyper_ycbcr_ref(:,:,2));
Cr_hyper_ref = double(hyper_ycbcr_ref(:,:,3));

delta_CbCr_cross = sqrt((Cb_film_ref - Cb_hyper_ref).^2 + ...
                        (Cr_film_ref - Cr_hyper_ref).^2);

sigma_small = 1;
sigma_large = 300;


% 5️⃣ DoG on CbCr film after msr and hyper no msr
DoG_Cb_film_ref  = imgaussfilt(Cb_film_ref, sigma_small) - imgaussfilt(Cb_film_ref, sigma_large);
DoG_Cr_film_ref  = imgaussfilt(Cr_film_ref, sigma_small) - imgaussfilt(Cr_film_ref, sigma_large);
DoG_Cb_hyper_ref = imgaussfilt(Cb_hyper_ref, sigma_small) - imgaussfilt(Cb_hyper_ref, sigma_large);
DoG_Cr_hyper_ref = imgaussfilt(Cr_hyper_ref, sigma_small) - imgaussfilt(Cr_hyper_ref, sigma_large);

DoG_diff_CbCr_cross = sqrt((DoG_Cb_film_ref - DoG_Cb_hyper_ref).^2 + ...
                           (DoG_Cr_film_ref - DoG_Cr_hyper_ref).^2);
DoG_diff_CbCr_cross = mat2gray(DoG_diff_CbCr_cross);

%  Visualization
figure('Name','Cross-Domain Comparisons (Film after MSR vs Hyper before MSR)');
subplot(1,1,1); imagesc(DoG_diff_CbCr_cross, [0 1]); axis image off; colormap(turbo); colorbar;
title('DoG on CbCr');

%%




%% ===========================================
% Advanced Preprocessing Pipeline for Film vs Hyperspectral Comparison
% Computes delta CbCr after each step
% ===========================================


%% Parameters
isGamma = true;
sigma_bilateral = 5;    % for bilateral filtering
sigma_large      = 300; % for DoG/edge extraction if needed

%

%% 5️⃣ Cross-Domain Linear Regression (MSR after vs Hyper before)
scales = [15,80,250];
film_msr = MSR(film_rgb, scales);

film_corr = zeros(size(film_rgb));
for c = 1:3
    x = reshape(film_msr(:,:,c), [], 1);
    y = reshape(hyper_rgb(:,:,c), [], 1);
    p = polyfit(x,y,1);
    film_corr(:,:,c) = reshape(polyval(p,x), size(film_msr(:,:,c)));
end

deltaCbCr_crossdomain = compute_deltaCbCr(film_corr, hyper_rgb);

% Visualization
figure('Name','ΔCbCr Comparison');
subplot(1,1,1); imagesc(deltaCbCr_crossdomain); axis image off; colormap(jet); colorbar; title('Cross-Domain Regression');





%%

%




%% -------------------------
%% Local function at the end of the script
function delta = compute_deltaCbCr(film, hyp)
    % Convert to YCbCr
    ycbcr_film = rgb2ycbcr(film);
    ycbcr_hyp  = rgb2ycbcr(hyp);

    % Extract Cb and Cr channels
    Cb_film = double(ycbcr_film(:,:,2));
    Cr_film = double(ycbcr_film(:,:,3));
    Cb_hyp  = double(ycbcr_hyp(:,:,2));
    Cr_hyp  = double(ycbcr_hyp(:,:,3));

    % Compute ΔCbCr
    delta = sqrt( (Cb_film - Cb_hyp).^2 + (Cr_film - Cr_hyp).^2 );
end


%%





