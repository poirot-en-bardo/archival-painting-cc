path_before = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat';
path_after  = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat';
% path_before = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat';
% path_after  = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat';

path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat';
% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_exp0.mat';
% path_film = '/Volumes/School/Thesis/data/captures/xyz_lab_rgb/film/cactus_led_fuji_underexp.mat';
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


%% --- Compute basic chroma maps (no MSR) ---


% ΔCbCr (YCbCr chroma)
film_ycbcr = rgb2ycbcr(film_rgb);
hyper_ycbcr = rgb2ycbcr(hyper_rgb);
delta_CbCr = sqrt((film_ycbcr(:,:,2) - hyper_ycbcr(:,:,2)).^2 + ...
                  (film_ycbcr(:,:,3) - hyper_ycbcr(:,:,3)).^2);


%% Compute several alternative w_cbcr (confidence for CbCr)
% Inputs assumed:
% Cb = film_ycbcr(:,:,2); Cr = film_ycbcr(:,:,3);
% use doubles
Cb = double(film_ycbcr(:,:,2));
Cr = double(film_ycbcr(:,:,3));

% parameters you can tune
patch = 2;        % local window for box / variance filters (odd recommended)
smooth_sigma = 2;  % smoothing of final weight map

% chroma magnitude (absolute chroma) — where C is strong, CbCr is informative
chroma_mag = sqrt(Cb.^2 + Cr.^2);
% w_cbcr_C = mat2gray(chroma_mag);
w_cbcr_C = mat2gray(imgaussfilt(chroma_mag, patch));
% w_cbcr_C = imgaussfilt(w_cbcr_C, smooth_sigma);


%% Visualize all candidate weights
figure('Name','Candidate w\_cbcr (confidence maps)','Position',[100 100 1200 700]);
subplot(1,1,1); imagesc(w_cbcr_C); axis image off; colorbar; colormap(turbo);title('C: chroma magnitude');

%%
film_ycbcr_nomsr  = rgb2ycbcr(film_rgb);
hyper_ycbcr_nomsr = rgb2ycbcr(hyper_rgb);

Cb_film_nomsr  = double(film_ycbcr_nomsr(:,:,2));
Cr_film_nomsr  = double(film_ycbcr_nomsr(:,:,3));

Cb_hyper_nomsr = double(hyper_ycbcr_nomsr(:,:,2));
Cr_hyper_nomsr = double(hyper_ycbcr_nomsr(:,:,3));

delta_CbCr_nomsr = sqrt((Cb_film_nomsr - Cb_hyper_nomsr).^2 + ...
                        (Cr_film_nomsr - Cr_hyper_nomsr).^2);

%%



%% --- Compute CbCr-based confidence map for material change ---

% Convert images to YCbCr
film_ycbcr  = rgb2ycbcr(film_rgb);
hyper_ycbcr = rgb2ycbcr(hyper_rgb);

% Extract Cb and Cr channels (double precision)
Cb_film   = double(film_ycbcr(:,:,2));
Cr_film   = double(film_ycbcr(:,:,3));
Cb_hyper  = double(hyper_ycbcr(:,:,2));
Cr_hyper  = double(hyper_ycbcr(:,:,3));

%% Parameters
patch = 2;         % Gaussian smoothing for local neighborhood
smooth_sigma = 2;  % optional extra smoothing

%% Compute chroma magnitude of the film (to weight important regions)
film_chroma_mag = sqrt(Cb_film.^2 + Cr_film.^2);

%% Compute difference in Cb/Cr between film and hyper image
delta_CbCr = sqrt((Cb_film - Cb_hyper).^2 + (Cr_film - Cr_hyper).^2);

%% Combine magnitude and difference
weighted_diff = film_chroma_mag .* delta_CbCr;

%% Smooth and normalize to [0,1]
w_cbcr_C = mat2gray(imgaussfilt(weighted_diff, patch));
w_cbcr_C = imgaussfilt(w_cbcr_C, smooth_sigma);

%% Visualize
figure('Name','CbCr Confidence Map (Material Change)','Position',[100 100 1200 700]);
imagesc(w_cbcr_C); axis image off; colorbar;
title('CbCr Confidence Map highlighting material change');
