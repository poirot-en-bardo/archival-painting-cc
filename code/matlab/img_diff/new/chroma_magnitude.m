ath_before = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat';
path_after  = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat';
% path_before = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat';
% path_after  = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat';

% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat';
% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_exp0.mat';
path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_underexp.mat';
% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/yoda_halogen_fuji_exp0.mat';
% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/yoda_halogen_fuji_overexp.mat';
% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/yoda_led_kodak_exp0.mat';

close all;

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

scales = [15, 80, 250];
film_msr = MSR(film_rgb, scales);
hyper_msr = MSR(hyper_rgb, scales);


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



%%
% Extract Lab channels before MSR
L_film   = Lab_film(:,:,1); 
a_film   = Lab_film(:,:,2); 
b_film   = Lab_film(:,:,3);

L_hyper  = Lab_after(:,:,1); 
a_hyper  = Lab_after(:,:,2); 
b_hyper  = Lab_after(:,:,3);



% Before MSR
ab_film_nomsr  = cat(3, a_film, b_film);
ab_hyper_nomsr = cat(3, a_hyper, b_hyper);

% After MSR
ab_film_msr  = cat(3, a_film_msr, b_film_msr);
ab_hyper_msr = cat(3, a_hyper_msr, b_hyper_msr);

%% --- Compute and visualize chroma magnitude-based w_cbcr (before & after MSR) ---

% --- Before MSR ---
film_ycbcr_nomsr  = rgb2ycbcr(film_rgb);
hyper_ycbcr_nomsr = rgb2ycbcr(hyper_rgb);

Cb_film_nomsr = double(film_ycbcr_nomsr(:,:,2));
Cr_film_nomsr = double(film_ycbcr_nomsr(:,:,3));
Cb_hyper_nomsr = double(hyper_ycbcr_nomsr(:,:,2));
Cr_hyper_nomsr = double(hyper_ycbcr_nomsr(:,:,3));

% Chroma magnitude maps (film + hyper)
chroma_mag_film_nomsr  = sqrt(Cb_film_nomsr.^2 + Cr_film_nomsr.^2);
chroma_mag_hyper_nomsr = sqrt(Cb_hyper_nomsr.^2 + Cr_hyper_nomsr.^2);

% --- After MSR ---
film_ycbcr_msr  = rgb2ycbcr(film_msr);
hyper_ycbcr_msr = rgb2ycbcr(hyper_msr);

Cb_film_msr = double(film_ycbcr_msr(:,:,2));
Cr_film_msr = double(film_ycbcr_msr(:,:,3));
Cb_hyper_msr = double(hyper_ycbcr_msr(:,:,2));
Cr_hyper_msr = double(hyper_ycbcr_msr(:,:,3));

chroma_mag_film_msr  = sqrt(Cb_film_msr.^2 + Cr_film_msr.^2);
chroma_mag_hyper_msr = sqrt(Cb_hyper_msr.^2 + Cr_hyper_msr.^2);

% --- Normalize for visualization ---
chroma_mag_film_nomsr  = mat2gray(chroma_mag_film_nomsr);
chroma_mag_hyper_nomsr = mat2gray(chroma_mag_hyper_nomsr);
chroma_mag_film_msr    = mat2gray(chroma_mag_film_msr);
chroma_mag_hyper_msr   = mat2gray(chroma_mag_hyper_msr);

% --- Visualization ---
figure('Name','Chroma Magnitude (w_{CbCr}) Before/After MSR','Position',[100 100 1200 700]);

subplot(2,2,1); imagesc(chroma_mag_film_nomsr); axis image off; colorbar; colormap(turbo);
title('Film chroma magnitude (before MSR)');

subplot(2,2,2); imagesc(chroma_mag_hyper_nomsr); axis image off; colorbar; colormap(turbo);
title('Hyper chroma magnitude (before MSR)');

subplot(2,2,3); imagesc(chroma_mag_film_msr); axis image off; colorbar; colormap(turbo);
title('Film chroma magnitude (after MSR)');

subplot(2,2,4); imagesc(chroma_mag_hyper_msr); axis image off; colorbar; colormap(turbo);
title('Hyper chroma magnitude (after MSR)');

%% --- Optional: compare difference in chroma magnitude ---
delta_chroma_mag_nomsr = abs(chroma_mag_film_nomsr - chroma_mag_hyper_nomsr);
delta_chroma_mag_msr   = abs(chroma_mag_film_msr   - chroma_mag_hyper_msr);

figure('Name','Δ Chroma Magnitude Before/After MSR','Position',[100 100 1000 400]);
subplot(1,2,1); imagesc(delta_chroma_mag_nomsr); axis image off; colorbar; colormap(turbo);
title('\Delta chroma magnitude (before MSR)');
subplot(1,2,2); imagesc(delta_chroma_mag_msr); axis image off; colorbar; colormap(turbo);
title('\Delta chroma magnitude (after MSR)');



%%
function YCbCr = prophoto2ycbcr_custom(RGB)
    % Convert ProPhoto RGB (D50, gamma encoded) → YCbCr (BT.709, D65)
    % RGB expected in [0,1]
    
    % 1) Convert to XYZ (D50) using your provided function
    XYZ_D50 = prophoto2xyz(RGB, true);

    % 2) Adapt D50 → D65 using Bradford transform
    M_Bradford = [ ...
        0.8951,  0.2664, -0.1614;
       -0.7502,  1.7135,  0.0367;
        0.0389, -0.0685,  1.0296];

    D50 = [0.9642; 1.0000; 0.8251];
    D65 = [0.95047; 1.00000; 1.08883];
    M_adapt = inv(M_Bradford) * diag(D65 ./ (M_Bradford * D50)) * M_Bradford;
    XYZ_D65 = reshape((reshape(XYZ_D50,[],3) * M_adapt.'), size(XYZ_D50));

    % 3) Convert XYZ(D65) → linear sRGB
    M_XYZ2sRGB = [ ...
         3.2406, -1.5372, -0.4986;
        -0.9689,  1.8758,  0.0415;
         0.0557, -0.2040,  1.0570];
    RGB_lin = reshape((reshape(XYZ_D65,[],3) * M_XYZ2sRGB.'), size(RGB));

    % Clip and gamma-encode (approx. sRGB)
    RGB_lin = max(min(RGB_lin,1),0);
    mask = (RGB_lin <= 0.0031308);
    RGB_gamma = zeros(size(RGB_lin));
    RGB_gamma(mask) = 12.92 * RGB_lin(mask);
    RGB_gamma(~mask) = 1.055 * RGB_lin(~mask).^(1/2.4) - 0.055;

    % 4) Apply standard BT.709 YCbCr conversion (full-range [0,1])
    T = [ ...
        0.2126,  0.7152,  0.0722;  % Y
       -0.1146, -0.3854,  0.5000;  % Cb
        0.5000, -0.4542, -0.0458]; % Cr

    YCbCr = reshape((reshape(RGB_gamma,[],3) * T.'), size(RGB));
    % Normalize to [0,1]
    YCbCr(:,:,2:3) = YCbCr(:,:,2:3) + 0.5;
    YCbCr = max(min(YCbCr,1),0);
end
