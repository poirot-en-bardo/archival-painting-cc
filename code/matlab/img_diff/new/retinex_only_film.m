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

%%

%% --- Setup ---
scales = [15, 80, 250]; % MSR scales
sigma_small = 1;
sigma_large = 750;
DoG = @(X) imgaussfilt(X, sigma_small) - imgaussfilt(X, sigma_large);

RGB_film  = film_data.RGB_lin_img;
RGB_hyper = painting_after.RGB_lin_img;

[h,w,~] = size(RGB_film);

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

%% --- Initialize storage ---
spaces = {'YCbCr','OpponentRGB'};
delta_all = cell(length(spaces),3);
DoG_all   = cell(length(spaces),3);

%% --- Standardized Opponent RGB (orthonormal) matrix ---
M_opp = [ 1/sqrt(2), -1/sqrt(2), 0;
          1/sqrt(6),  1/sqrt(6), -2/sqrt(6);
          1/sqrt(3),  1/sqrt(3),  1/sqrt(3) ];  % De Valois & De Valois, 1988

%% --- Loop over MSR cases ---
for k = 1:3
    rgbF = max(cases{k}{2},0); % ensure positivity
    rgbH = max(cases{k}{3},0);

    %% --- 1) YCbCr ---
    YCbCr_F = cat(3, ...
        0.299*rgbF(:,:,1) + 0.587*rgbF(:,:,2) + 0.114*rgbF(:,:,3), ...
        0.5*(rgbF(:,:,3) - (0.299*rgbF(:,:,1) + 0.587*rgbF(:,:,2) + 0.114*rgbF(:,:,3)))./(1-0.114+eps), ...
        0.5*(rgbF(:,:,1) - (0.299*rgbF(:,:,1) + 0.587*rgbF(:,:,2) + 0.114*rgbF(:,:,3)))./(1-0.299+eps));

    YCbCr_H = cat(3, ...
        0.299*rgbH(:,:,1) + 0.587*rgbH(:,:,2) + 0.114*rgbH(:,:,3), ...
        0.5*(rgbH(:,:,3) - (0.299*rgbH(:,:,1) + 0.587*rgbH(:,:,2) + 0.114*rgbH(:,:,3)))./(1-0.114+eps), ...
        0.5*(rgbH(:,:,1) - (0.299*rgbH(:,:,1) + 0.587*rgbH(:,:,2) + 0.114*rgbH(:,:,3)))./(1-0.299+eps));

    delta_all{1,k} = mat2gray(sqrt((YCbCr_F(:,:,2)-YCbCr_H(:,:,2)).^2 + ...
                                    (YCbCr_F(:,:,3)-YCbCr_H(:,:,3)).^2));
    DoG_all{1,k}   = mat2gray(sqrt((DoG(YCbCr_F(:,:,2))-DoG(YCbCr_H(:,:,2))).^2 + ...
                                    (DoG(YCbCr_F(:,:,3))-DoG(YCbCr_H(:,:,3))).^2));

    %% --- 2) Standardized Opponent RGB ---
    Opp_F = reshape((reshape(rgbF,[],3) * M_opp.'), h, w, 3);
    Opp_H = reshape((reshape(rgbH,[],3) * M_opp.'), h, w, 3);

    % Compute differences using chromatic channels only (O1, O2)
    delta_all{2,k} = mat2gray(sqrt(sum((Opp_F(:,:,1:2)-Opp_H(:,:,1:2)).^2,3)));
    DoG_all{2,k}   = mat2gray(sqrt(sum((DoG(Opp_F(:,:,1:2))-DoG(Opp_H(:,:,1:2))).^2,3)));
end

%% --- Visualization ---
for s = 1:length(spaces)
    figure('Name',spaces{s},'Position',[50 50 1400 700]);
    for k = 1:3
        subplot(2,3,k); imagesc(delta_all{s,k}); axis image off; colorbar; colormap(jet);
        title([spaces{s} ' Δ – ' cases{k}{1}]);
        subplot(2,3,3+k); imagesc(DoG_all{s,k}); axis image off; colorbar; colormap(turbo);
        title([spaces{s} ' DoG – ' cases{k}{1}]);
    end
end


















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

function out = MSRCR(img, scales, alpha, beta, gain, offset)
% MSRCR  Multi-Scale Retinex with Color Restoration
%
%   out = MSRCR(img, scales, alpha, beta, gain, offset)
%
%   img    : HxWx3 RGB image (double, recommended range 0–1)
%   scales : array of Gaussian sigmas, e.g. [15 80 250]
%   alpha  : MSRCR parameter (default 125)
%   beta   : MSRCR parameter (default 46)
%   gain   : output global gain (default 1)
%   offset : output offset (default 0)
%
%   Reference: Rahman, Jobson, Woodell — MSRCR (1996, 1997, 2004)

    if nargin < 3, alpha  = 125; end
    if nargin < 4, beta   = 46;  end
    if nargin < 5, gain   = 1;   end
    if nargin < 6, offset = 0;   end

    img = double(img);
    img(img < 1e-6) = 1e-6;     % avoid log(0)

    % --- Multi-scale Retinex ---
    N = numel(scales);
    msr = zeros(size(img));

    for c = 1:3
        I = img(:,:,c);
        ret = zeros(size(I));

        for s = scales
            blurred = imgaussfilt(I, s);
            ret = ret + (log(I) - log(blurred));
        end

        msr(:,:,c) = ret / N;
    end

    % --- Color Restoration Function (CRF) ---
    %    C = beta * ( log( alpha * I ) - log( R+G+B ) )
    sumRGB = img(:,:,1) + img(:,:,2) + img(:,:,3);
    sumRGB(sumRGB < 1e-6) = 1e-6;

    CRF = beta * ( log(alpha * img) - log(sumRGB) );

    % --- MSRCR ---
    out = gain * (CRF .* msr) + offset;

    % Optional global normalization to [0,1]
    out = out - min(out(:));
    out = out / max(out(:) + eps);

end


scales = [1, 20, 200];
film_msr = MSR(film_rgb, scales);
hyper_msr = MSR(hyper_rgb, scales);
% hyper_msr = hyper_rgb;

% Intermediate visualization
% figure;
% subplot(2,2,1); imshow(film_rgb); title('Original Film');
% subplot(2,2,2); imshow(film_msr); title('Film after MSR');
% subplot(2,2,3); imshow(hyper_rgb); title('Original Hyper');
% subplot(2,2,4); imshow(hyper_msr); title('Hyper after MSR');



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
% figure;
% subplot(2,3,1); imshow(L_film_msr,[]); title('Film L (after MSR)');
% subplot(2,3,2); imshow(a_film_msr,[]); title('Film a (after MSR)');
% subplot(2,3,3); imshow(b_film_msr,[]); title('Film b (after MSR)');
% subplot(2,3,4); imshow(L_hyper_msr,[]); title('Hyper L (after MSR)');
% subplot(2,3,5); imshow(a_hyper_msr,[]); title('Hyper a (after MSR)');
% subplot(2,3,6); imshow(b_hyper_msr,[]); title('Hyper b (after MSR)');

%%
% Extract Lab channels before MSR
L_film   = Lab_film(:,:,1); 
a_film   = Lab_film(:,:,2); 
b_film   = Lab_film(:,:,3);

L_hyper  = Lab_after(:,:,1); 
a_hyper  = Lab_after(:,:,2); 
b_hyper  = Lab_after(:,:,3);

% % Visualization of Lab channels before MSR
% figure;
% subplot(2,3,1); imshow(L_film, []); title('Film L (before MSR)');
% subplot(2,3,2); imshow(a_film, []); title('Film a (before MSR)');
% subplot(2,3,3); imshow(b_film, []); title('Film b (before MSR)');
% subplot(2,3,4); imshow(L_hyper, []); title('Hyper L (before MSR)');
% subplot(2,3,5); imshow(a_hyper, []); title('Hyper a (before MSR)');
% subplot(2,3,6); imshow(b_hyper, []); title('Hy per b (before MSR)');

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




%% Compute delta_CbCr (before and after MSR)

% --- Before MSR ---
% film_ycbcr_nomsr  = rgb2ycbcr(film_rgb);
% hyper_ycbcr_nomsr = rgb2ycbcr(hyper_rgb);
% 
% Cb_film_nomsr  = double(film_ycbcr_nomsr(:,:,2));
% Cr_film_nomsr  = double(film_ycbcr_nomsr(:,:,3));
% 
% Cb_hyper_nomsr = double(hyper_ycbcr_nomsr(:,:,2));
% Cr_hyper_nomsr = double(hyper_ycbcr_nomsr(:,:,3));
% 
% delta_CbCr_nomsr = sqrt((Cb_film_nomsr - Cb_hyper_nomsr).^2 + ...
%                         (Cr_film_nomsr - Cr_hyper_nomsr).^2);
% 
% % --- After MSR ---
% film_ycbcr_msr  = rgb2ycbcr(film_msr);
% hyper_ycbcr_msr = rgb2ycbcr(hyper_rgb); %change
% 
% Cb_film_msr  = double(film_ycbcr_msr(:,:,2));
% Cr_film_msr  = double(film_ycbcr_msr(:,:,3));
% 
% Cb_hyper_msr = double(hyper_ycbcr_msr(:,:,2));
% Cr_hyper_msr = double(hyper_ycbcr_msr(:,:,3));
% 
% delta_CbCr_msr = sqrt((Cb_film_msr - Cb_hyper_msr).^2 + ...
%                       (Cr_film_msr - Cr_hyper_msr).^2);
% 
% % Visualization
% figure;
% subplot(1,2,1);
% imagesc(delta_CbCr_nomsr); axis image off; colorbar;
% colormap(jet); title('\Delta CbCr (before MSR)');
% 
% subplot(1,2,2);
% imagesc(delta_CbCr_msr); axis image off; colorbar;
% colormap(jet); title('\Delta CbCr (after MSR)');

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
% subplot(2,3,1); imshow(Y_film_nomsr,[]); title('Film Y (before MSR)');
% subplot(2,3,2); imshow(Cb_film_nomsr,[]); title('Film Cb (before MSR)');
% subplot(2,3,3); imshow(Cr_film_nomsr,[]); title('Film Cr (before MSR)');
% subplot(2,3,4); imshow(Y_hyper_nomsr,[]); title('Hyper Y (before MSR)');
% subplot(2,3,5); imshow(Cb_hyper_nomsr,[]); title('Hyper Cb (before MSR)');
% subplot(2,3,6); imshow(Cr_hyper_nomsr,[]); title('Hyper Cr (before MSR)');

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
% subplot(2,3,1); imshow(Y_film_msr,[]); title('Film Y (after MSR)');
% subplot(2,3,2); imshow(Cb_film_msr,[]); title('Film Cb (after MSR)');
% subplot(2,3,3); imshow(Cr_film_msr,[]); title('Film Cr (after MSR)');
% subplot(2,3,4); imshow(Y_hyper_msr,[]); title('Hyper Y (after MSR)');
% subplot(2,3,5); imshow(Cb_hyper_msr,[]); title('Hyper Cb (after MSR)');
% subplot(2,3,6); imshow(Cr_hyper_msr,[]); title('Hyper Cr (after MSR)');


% %%
% %% DoG on CbCr channels (before and after MSR)
% 
% sigma_small = 1;
% sigma_large = 300;
% 
% % --- Before MSR ---
% DoG_Cb_film_nomsr  = imgaussfilt(Cb_film_nomsr, sigma_small) - imgaussfilt(Cb_film_nomsr, sigma_large);
% DoG_Cr_film_nomsr  = imgaussfilt(Cr_film_nomsr, sigma_small) - imgaussfilt(Cr_film_nomsr, sigma_large);
% DoG_Cb_hyper_nomsr = imgaussfilt(Cb_hyper_nomsr, sigma_small) - imgaussfilt(Cb_hyper_nomsr, sigma_large);
% DoG_Cr_hyper_nomsr = imgaussfilt(Cr_hyper_nomsr, sigma_small) - imgaussfilt(Cr_hyper_nomsr, sigma_large);
% 
% DoG_diff_CbCr_nomsr = sqrt( (DoG_Cb_film_nomsr - DoG_Cb_hyper_nomsr).^2 + ...
%                             (DoG_Cr_film_nomsr - DoG_Cr_hyper_nomsr).^2 );
% 
% % --- After MSR ---
% DoG_Cb_film_msr  = imgaussfilt(Cb_film_msr, sigma_small) - imgaussfilt(Cb_film_msr, sigma_large);
% DoG_Cr_film_msr  = imgaussfilt(Cr_film_msr, sigma_small) - imgaussfilt(Cr_film_msr, sigma_large);
% DoG_Cb_hyper_msr = imgaussfilt(Cb_hyper_msr, sigma_small) - imgaussfilt(Cb_hyper_msr, sigma_large);
% DoG_Cr_hyper_msr = imgaussfilt(Cr_hyper_msr, sigma_small) - imgaussfilt(Cr_hyper_msr, sigma_large);
% 
% DoG_diff_CbCr_msr = sqrt( (DoG_Cb_film_msr - DoG_Cb_hyper_msr).^2 + ...
%                           (DoG_Cr_film_msr - DoG_Cr_hyper_msr).^2 );
% 
% % Normalize for visualization
% DoG_diff_CbCr_nomsr = mat2gray(DoG_diff_CbCr_nomsr);
% DoG_diff_CbCr_msr   = mat2gray(DoG_diff_CbCr_msr);
% 
% % --- Visualization ---
% figure;
% subplot(1,2,1);
% imagesc(DoG_diff_CbCr_nomsr, [0 1]); axis image off; colormap(turbo); colorbar;
% title('DoG on CbCr (before MSR)');
% 
% subplot(1,2,2);
% imagesc(DoG_diff_CbCr_msr, [0 1]); axis image off; colormap(turbo); colorbar;
% title('DoG on CbCr (after MSR)');

%%







%% 🔁 Cross comparison: Film after MSR vs Hyper before MSR
% 
% % Define your "reference" pair
% film_ref   = film_msr;     % corrected film
% hyper_ref  = hyper_rgb;    % uncorrected hyperspectral RGB
% 
% % Convert both to Lab
% isGamma = true;
% XYZ_film_ref  = prophoto2xyz(film_ref, isGamma);
% XYZ_hyper_ref = prophoto2xyz(hyper_ref, isGamma);
% 
% film_lab_ref  = xyz2lab_custom(reshape(XYZ_film_ref,[],3));
% hyper_lab_ref = xyz2lab_custom(reshape(XYZ_hyper_ref,[],3));
% 
% film_lab_ref  = reshape(film_lab_ref, size(film_ref));
% hyper_lab_ref = reshape(hyper_lab_ref, size(hyper_ref));
% 
% % Extract channels
% a_film_ref = film_lab_ref(:,:,2);
% b_film_ref = film_lab_ref(:,:,3);
% a_hyper_ref = hyper_lab_ref(:,:,2);
% b_hyper_ref = hyper_lab_ref(:,:,3);
% 
% % 3️⃣ ΔCbCr
% film_ycbcr_ref  = rgb2ycbcr(film_ref);
% hyper_ycbcr_ref = rgb2ycbcr(hyper_ref);
% 
% Cb_film_ref = double(film_ycbcr_ref(:,:,2));
% Cr_film_ref = double(film_ycbcr_ref(:,:,3));
% Cb_hyper_ref = double(hyper_ycbcr_ref(:,:,2));
% Cr_hyper_ref = double(hyper_ycbcr_ref(:,:,3));
% 
% delta_CbCr_cross = sqrt((Cb_film_ref - Cb_hyper_ref).^2 + ...
%                         (Cr_film_ref - Cr_hyper_ref).^2);
% 
% sigma_small = 1;
% sigma_large = 300;
% 
% 
% % 5️⃣ DoG on CbCr film after msr and hyper no msr
% DoG_Cb_film_ref  = imgaussfilt(Cb_film_ref, sigma_small) - imgaussfilt(Cb_film_ref, sigma_large);
% DoG_Cr_film_ref  = imgaussfilt(Cr_film_ref, sigma_small) - imgaussfilt(Cr_film_ref, sigma_large);
% DoG_Cb_hyper_ref = imgaussfilt(Cb_hyper_ref, sigma_small) - imgaussfilt(Cb_hyper_ref, sigma_large);
% DoG_Cr_hyper_ref = imgaussfilt(Cr_hyper_ref, sigma_small) - imgaussfilt(Cr_hyper_ref, sigma_large);
% 
% DoG_diff_CbCr_cross = sqrt((DoG_Cb_film_ref - DoG_Cb_hyper_ref).^2 + ...
%                            (DoG_Cr_film_ref - DoG_Cr_hyper_ref).^2);
% DoG_diff_CbCr_cross = mat2gray(DoG_diff_CbCr_cross);
% 
% %  Visualization
% figure('Name','Cross-Domain Comparisons (Film after MSR vs Hyper before MSR)');
% subplot(1,1,1); imagesc(DoG_diff_CbCr_cross, [0 1]); axis image off; colormap(turbo); colorbar;
% title('DoG on CbCr');
% 





%% --- ALL ---


% Scales for MSR
scales = [15, 80, 250];

% Apply MSR only to FILM
film_msr = MSR(film_rgb, scales);

% Hyper stays uncorrected
hyper_msr = hyper_rgb;

% %% --- Visualize MSR effect ---
% figure;
% subplot(1,3,1); imshow(film_rgb); title('Film (original, linear)');
% subplot(1,3,2); imshow(film_msr); title('Film (after MSR)');
% subplot(1,3,3); imshow(hyper_rgb); title('Hyper (uncorrected)');
% 
% 
% %% --- Step 2: Convert both to XYZ, then Lab ---
% isGamma = false; % linear RGB → XYZ (no gamma correction)
% 
% XYZ_film_msr  = prophoto2xyz(film_msr, isGamma);
% XYZ_hyper_msr = prophoto2xyz(hyper_msr, isGamma);
% 
% film_lab_msr  = xyz2lab_custom(reshape(XYZ_film_msr, [], 3));
% hyper_lab_msr = xyz2lab_custom(reshape(XYZ_hyper_msr, [], 3));
% 
% film_lab_msr  = reshape(film_lab_msr, size(film_msr));
% hyper_lab_msr = reshape(hyper_lab_msr, size(hyper_msr));
% 
% 
% %% --- Step 1: Apply MSR to linear film RGB ---
% scales = [15, 80, 250];
% film_msr = MSR(film_rgb, scales);   % linear RGB
% 
% %% --- Step 2: Convert to YCbCr ---
% % Use MATLAB rgb2ycbcr on linear images
% film_ycbcr_nomsr  = rgb2ycbcr(film_rgb);
% film_ycbcr_msr    = rgb2ycbcr(film_msr);
% hyper_ycbcr       = rgb2ycbcr(hyper_rgb);
% 
% %% --- Step 3: ΔCbCr computations ---
% % No MSR (film original vs hyper)
% delta_CbCr_nomsr = sqrt((double(film_ycbcr_nomsr(:,:,2)) - double(hyper_ycbcr(:,:,2))).^2 + ...
%                         (double(film_ycbcr_nomsr(:,:,3)) - double(hyper_ycbcr(:,:,3))).^2);
% 
% % MSR on film only
% delta_CbCr_msr = sqrt((double(film_ycbcr_msr(:,:,2)) - double(hyper_ycbcr(:,:,2))).^2 + ...
%                       (double(film_ycbcr_msr(:,:,3)) - double(hyper_ycbcr(:,:,3))).^2);
% 
% %% --- Step 4: DoG on CbCr ---
% sigma_small = 1;
% sigma_large = 300;
% 
% % Helper for DoG
% DoG = @(X) imgaussfilt(X, sigma_small) - imgaussfilt(X, sigma_large);
% 
% % No MSR
% DoG_diff_CbCr_nomsr = sqrt( (DoG(double(film_ycbcr_nomsr(:,:,2))) - DoG(double(hyper_ycbcr(:,:,2)))).^2 + ...
%                             (DoG(double(film_ycbcr_nomsr(:,:,3))) - DoG(double(hyper_ycbcr(:,:,3)))).^2 );
% 
% % MSR on film only
% DoG_diff_CbCr_msr = sqrt( (DoG(double(film_ycbcr_msr(:,:,2))) - DoG(double(hyper_ycbcr(:,:,2)))).^2 + ...
%                           (DoG(double(film_ycbcr_msr(:,:,3))) - DoG(double(hyper_ycbcr(:,:,3)))).^2 );
% 
% % Normalize for display
% DoG_diff_CbCr_nomsr = mat2gray(DoG_diff_CbCr_nomsr);
% DoG_diff_CbCr_msr   = mat2gray(DoG_diff_CbCr_msr);
% 
% % --- Step 5: Visualization ---
% figure('Name','ΔCbCr and DoG CbCr Comparisons');
% 
% subplot(2,2,1);
% imagesc(delta_CbCr_nomsr); axis image off; colorbar; colormap(jet);
% title('\Delta CbCr (before MSR)');
% 
% subplot(2,2,2);
% imagesc(delta_CbCr_msr); axis image off; colorbar; colormap(jet);
% title('\Delta CbCr (after MSR on film)');
% 
% subplot(2,2,3);
% imagesc(DoG_diff_CbCr_nomsr); axis image off; colorbar; colormap(turbo);
% title('DoG CbCr (before MSR)');
% 
% subplot(2,2,4);
% imagesc(DoG_diff_CbCr_msr); axis image off; colorbar; colormap(turbo);
% title('DoG CbCr (after MSR on film)');
% 





%%






%% --- Step 1: Apply MSR ---
scales = [100, 700, 1050];
% scales = [15, 80, 250]; % MSR scales


film_msr  = MSR(film_rgb, scales);
hyper_msr = MSR(hyper_rgb, scales);

% --- Step 2: Convert to YCbCr (linear RGB) ---
film_ycbcr_nomsr  = rgb2ycbcr(film_rgb);
film_ycbcr_msr    = rgb2ycbcr(film_msr);
hyper_ycbcr       = rgb2ycbcr(hyper_rgb);
hyper_ycbcr_msr   = rgb2ycbcr(hyper_msr);

% --- Step 3: ΔCbCr Computations ---
% Case 1: No MSR
delta_CbCr_nomsr = sqrt((double(film_ycbcr_nomsr(:,:,2)) - double(hyper_ycbcr(:,:,2))).^2 + ...
                        (double(film_ycbcr_nomsr(:,:,3)) - double(hyper_ycbcr(:,:,3))).^2);

% Case 2: Both MSR
delta_CbCr_both = sqrt((double(film_ycbcr_msr(:,:,2)) - double(hyper_ycbcr_msr(:,:,2))).^2 + ...
                       (double(film_ycbcr_msr(:,:,3)) - double(hyper_ycbcr_msr(:,:,3))).^2);

% Case 3: Film-only MSR (cross-domain)
delta_CbCr_cross = sqrt((double(film_ycbcr_msr(:,:,2)) - double(hyper_ycbcr(:,:,2))).^2 + ...
                        (double(film_ycbcr_msr(:,:,3)) - double(hyper_ycbcr(:,:,3))).^2);
%%
% --- Step 4: DoG on CbCr ---
sigma_small = 1;
sigma_large = 500;
DoG = @(X) imgaussfilt(X, sigma_small) - imgaussfilt(X, sigma_large);

% Case 1: No MSR
DoG_diff_CbCr_nomsr = sqrt((DoG(double(film_ycbcr_nomsr(:,:,2))) - DoG(double(hyper_ycbcr(:,:,2)))).^2 + ...
                           (DoG(double(film_ycbcr_nomsr(:,:,3))) - DoG(double(hyper_ycbcr(:,:,3)))).^2);

% Case 2: Both MSR
DoG_diff_CbCr_both = sqrt((DoG(double(film_ycbcr_msr(:,:,2))) - DoG(double(hyper_ycbcr_msr(:,:,2)))).^2 + ...
                          (DoG(double(film_ycbcr_msr(:,:,3))) - DoG(double(hyper_ycbcr_msr(:,:,3)))).^2);

% Case 3: Film-only MSR
DoG_diff_CbCr_cross = sqrt((DoG(double(film_ycbcr_msr(:,:,2))) - DoG(double(hyper_ycbcr(:,:,2)))).^2 + ...
                           (DoG(double(film_ycbcr_msr(:,:,3))) - DoG(double(hyper_ycbcr(:,:,3)))).^2);

% Normalize for display
DoG_diff_CbCr_nomsr = mat2gray(DoG_diff_CbCr_nomsr);
DoG_diff_CbCr_both  = mat2gray(DoG_diff_CbCr_both);
DoG_diff_CbCr_cross = mat2gray(DoG_diff_CbCr_cross);

%%
% --- Step 5: Visualization ---
figure('Name','ΔCbCr and DoG CbCr Comparisons (3 MSR Cases)', 'Position',[100 100 1400 800]);

subplot(2,3,1);
imagesc(delta_CbCr_nomsr); axis image off; colorbar; colormap(parula); 
title('\Delta CbCr — No MSR');

subplot(2,3,2);
imagesc(delta_CbCr_both); axis image off; colorbar; colormap(parula); 
title('\Delta CbCr — Both MSR');

subplot(2,3,3);
imagesc(delta_CbCr_cross); axis image off; colorbar; colormap(parula); 
title('\Delta CbCr — Film-only MSR');

subplot(2,3,4);
imagesc(DoG_diff_CbCr_nomsr); axis image off; colorbar; colormap(parula); 
% clim([0, 0.5]);
title('DoG CbCr — No MSR');

subplot(2,3,5);
imagesc(DoG_diff_CbCr_both); axis image off; colorbar; colormap(turbo); 
% clim([0, 0.5]);
title('DoG CbCr — Both MSR');

subplot(2,3,6);
imagesc(DoG_diff_CbCr_cross); axis image off; colorbar; colormap(jet); 
% clim([0, 0.5]);
title('DoG CbCr — Film-only MSR');

sgtitle('ΔCbCr and DoG CbCr across MSR Conditions');










%%







%% --- ΔCbCr and DoG CbCr for Painting Before vs After ---

% Extract RGB images
painting_before_rgb = painting_before.RGB_img;
painting_after_rgb  = painting_after.RGB_img;

% Convert both to YCbCr
ycbcr_before = rgb2ycbcr(painting_before_rgb);
ycbcr_after  = rgb2ycbcr(painting_after_rgb);

Cb_before = double(ycbcr_before(:,:,2));
Cr_before = double(ycbcr_before(:,:,3));

Cb_after  = double(ycbcr_after(:,:,2));
Cr_after  = double(ycbcr_after(:,:,3));

% ΔCbCr
deltaCbCr_painting = sqrt((Cb_after - Cb_before).^2 + (Cr_after - Cr_before).^2);

% DoG on CbCr channels
sigma_small = 1;
sigma_large = 750;
DoG = @(X) imgaussfilt(X, sigma_small) - imgaussfilt(X, sigma_large);

DoG_diffCbCr_painting = sqrt((DoG(Cb_after) - DoG(Cb_before)).^2 + ...
                             (DoG(Cr_after) - DoG(Cr_before)).^2);

% Normalize for display
deltaCbCr_painting    = mat2gray(deltaCbCr_painting);
DoG_diffCbCr_painting = mat2gray(DoG_diffCbCr_painting);

%%
% --- Visualization ---
figure('Name','ΔCbCr and DoG CbCr — Painting Before vs After');
subplot(1,2,1);
imagesc(deltaCbCr_painting); axis image off; colorbar; colormap(parula);
caxis([0 0.5]);
title('ΔCbCr (Painting Before vs After)');

subplot(1,2,2);
imagesc(DoG_diffCbCr_painting); axis image off; colorbar; colormap(parula);
caxis([0 0.5]);
title('DoG CbCr (Painting Before vs After)');






%%
%% --- Painting ΔCbCr + DoG with Chromatic Adaptation (D50→D65) ---

% Extract XYZ images (already linear, D50)
XYZ_before = painting_before.XYZ_img;
XYZ_after  = painting_after.XYZ_img;

%% --- Chromatic adaptation D50 -> D65 ---
M_Bradford = [ 0.8951  0.2664 -0.1614;
              -0.7502  1.7135  0.0367;
               0.0389 -0.0685  1.0296 ];

D50_white = [0.9642; 1.0000; 0.8251];
D65_white = [0.9504; 1.0000; 1.0889];

M_adapt = inv(M_Bradford) * diag(M_Bradford*D65_white ./ (M_Bradford*D50_white)) * M_Bradford;

[h,w,~] = size(XYZ_before);
XYZ_before_D65 = reshape((M_adapt * reshape(XYZ_before,[],3)')', h,w,3);
XYZ_after_D65  = reshape((M_adapt * reshape(XYZ_after,[],3)')', h,w,3);

%% --- Convert adapted XYZ → linear sRGB ---
M_xyz2srgb = [ 3.2406 -1.5372 -0.4986;
              -0.9689  1.8758  0.0415;
               0.0557 -0.2040  1.0570 ];

rgb_before_lin = reshape((M_xyz2srgb * reshape(XYZ_before_D65,[],3)')', h,w,3);
rgb_after_lin  = reshape((M_xyz2srgb * reshape(XYZ_after_D65,[],3)')', h,w,3);

rgb_before_lin = max(rgb_before_lin,0);
rgb_after_lin  = max(rgb_after_lin,0);

%% --- Convert linear RGB → YCbCr ---
ycbcr_before = rgb2ycbcr(rgb_before_lin);
ycbcr_after  = rgb2ycbcr(rgb_after_lin);

Cb_before = double(ycbcr_before(:,:,2));
Cr_before = double(ycbcr_before(:,:,3));
Cb_after  = double(ycbcr_after(:,:,2));
Cr_after  = double(ycbcr_after(:,:,3));

%% --- ΔCbCr ---
deltaCbCr = sqrt((Cb_after - Cb_before).^2 + (Cr_after - Cr_before).^2);

%% --- DoG on CbCr ---
sigma_small = 1;
sigma_large = 750;
DoG = @(X) imgaussfilt(X,sigma_small) - imgaussfilt(X,sigma_large);

DoG_diffCbCr = sqrt((DoG(Cb_after) - DoG(Cb_before)).^2 + ...
                    (DoG(Cr_after) - DoG(Cr_before)).^2);

%% --- Normalize ---
deltaCbCr    = mat2gray(deltaCbCr);
DoG_diffCbCr = mat2gray(DoG_diffCbCr);

%% --- Visualization ---
figure('Name','Painting ΔCbCr and DoG CbCr (D50→D65 adapted)');
subplot(1,2,1);
imagesc(deltaCbCr); axis image off; colorbar; colormap(parula); caxis([0 0.5]);
title('ΔCbCr (Before vs After)');

subplot(1,2,2);
imagesc(DoG_diffCbCr); axis image off; colorbar; colormap(parula); caxis([0 0.5]);
title('DoG CbCr (Before vs After)');











%% retinex linear rgb


scales = [15, 80, 250]; % MSR scales
sigma_small = 1;
sigma_large = 300;
DoG = @(X) imgaussfilt(X, sigma_small) - imgaussfilt(X, sigma_large);

RGB_film = film_data.RGB_lin_img;
RGB_hyper = painting_after.RGB_lin_img;

[h,w,~] = size(RGB_film);

%% --- Prepare 3 MSR cases ---
film_case1  = RGB_film;             % no MSR
hyper_case1 = RGB_hyper;

film_case2  = MSR(RGB_film, scales); % both MSR
hyper_case2 = MSR(RGB_hyper, scales);

film_case3  = MSR(RGB_film, scales); % film only
hyper_case3 = RGB_hyper;

cases = {
    {'No MSR',        film_case1, hyper_case1};
    {'Both MSR',      film_case2, hyper_case2};
    {'Film only MSR', film_case3, hyper_case3};
};

%% --- Convert RGB → linear sRGB → YCbCr, compute ΔCbCr + DoG
deltaCbCr_all = cell(3,1);
DoGCbCr_all   = cell(3,1);

M_rgb2srgb = [ 3.2406, -1.5372, -0.4986;
              -0.9689,  1.8758,  0.0415;
               0.0557, -0.2040,  1.0570 ];

linear_rgb_to_ycbcr = @(rgb) cat(3, ...
    0.299*rgb(:,:,1) + 0.587*rgb(:,:,2) + 0.114*rgb(:,:,3), ...
    0.5*(rgb(:,:,3) - (0.299*rgb(:,:,1) + 0.587*rgb(:,:,2) + 0.114*rgb(:,:,3)))./(1-0.114+eps), ...
    0.5*(rgb(:,:,1) - (0.299*rgb(:,:,1) + 0.587*rgb(:,:,2) + 0.114*rgb(:,:,3)))./(1-0.299+eps) ...
);

for k = 1:3
    rgbF = cases{k}{2};
    rgbH = cases{k}{3};

    % Ensure positive values for sRGB conversion
    rgbF_lin = max(rgbF,0);
    rgbH_lin = max(rgbH,0);

    % Compute YCbCr
    YCbCrF = linear_rgb_to_ycbcr(rgbF_lin);
    YCbCrH = linear_rgb_to_ycbcr(rgbH_lin);

    % ΔCbCr
    deltaCbCr_all{k} = mat2gray(sqrt((YCbCrF(:,:,2)-YCbCrH(:,:,2)).^2 + ...
                                     (YCbCrF(:,:,3)-YCbCrH(:,:,3)).^2));

    % DoG on CbCr
    DoGCbCr_all{k} = mat2gray(sqrt((DoG(YCbCrF(:,:,2))-DoG(YCbCrH(:,:,2))).^2 + ...
                                    (DoG(YCbCrF(:,:,3))-DoG(YCbCrH(:,:,3))).^2));
end

%% --- Visualization
figure('Name','ΔCbCr and DoG CbCr (3 MSR cases)','Position',[100 100 1200 500]);
for k = 1:3
    subplot(2,3,k); imagesc(deltaCbCr_all{k}); axis image off; colorbar; colormap(jet);
    title(['ΔCbCr – ' cases{k}{1}]);
    subplot(2,3,3+k); imagesc(DoGCbCr_all{k}); axis image off; colorbar; colormap(turbo);
    title(['DoG CbCr – ' cases{k}{1}]);
end














%%
%all


%% --- Setup ---
scales = [15, 80, 250]; % MSR scales
sigma_small = 1;
sigma_large = 300;
DoG = @(X) imgaussfilt(X, sigma_small) - imgaussfilt(X, sigma_large);

RGB_film  = film_data.RGB_lin_img;
RGB_hyper = painting_after.RGB_lin_img;

[h,w,~] = size(RGB_film);

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

%% --- Initialize storage ---
spaces = {'YCbCr','OpponentRGB'};
delta_all = cell(length(spaces),3);
DoG_all   = cell(length(spaces),3);

%% --- Standardized Opponent RGB (orthonormal) matrix ---
M_opp = [ 1/sqrt(2), -1/sqrt(2), 0;
          1/sqrt(6),  1/sqrt(6), -2/sqrt(6);
          1/sqrt(3),  1/sqrt(3),  1/sqrt(3) ];  % De Valois & De Valois, 1988

%% --- Loop over MSR cases ---
for k = 1:3
    rgbF = max(cases{k}{2},0); % ensure positivity
    rgbH = max(cases{k}{3},0);

    %% --- 1) YCbCr ---
    YCbCr_F = cat(3, ...
        0.299*rgbF(:,:,1) + 0.587*rgbF(:,:,2) + 0.114*rgbF(:,:,3), ...
        0.5*(rgbF(:,:,3) - (0.299*rgbF(:,:,1) + 0.587*rgbF(:,:,2) + 0.114*rgbF(:,:,3)))./(1-0.114+eps), ...
        0.5*(rgbF(:,:,1) - (0.299*rgbF(:,:,1) + 0.587*rgbF(:,:,2) + 0.114*rgbF(:,:,3)))./(1-0.299+eps));

    YCbCr_H = cat(3, ...
        0.299*rgbH(:,:,1) + 0.587*rgbH(:,:,2) + 0.114*rgbH(:,:,3), ...
        0.5*(rgbH(:,:,3) - (0.299*rgbH(:,:,1) + 0.587*rgbH(:,:,2) + 0.114*rgbH(:,:,3)))./(1-0.114+eps), ...
        0.5*(rgbH(:,:,1) - (0.299*rgbH(:,:,1) + 0.587*rgbH(:,:,2) + 0.114*rgbH(:,:,3)))./(1-0.299+eps));

    delta_all{1,k} = mat2gray(sqrt((YCbCr_F(:,:,2)-YCbCr_H(:,:,2)).^2 + ...
                                    (YCbCr_F(:,:,3)-YCbCr_H(:,:,3)).^2));
    DoG_all{1,k}   = mat2gray(sqrt((DoG(YCbCr_F(:,:,2))-DoG(YCbCr_H(:,:,2))).^2 + ...
                                    (DoG(YCbCr_F(:,:,3))-DoG(YCbCr_H(:,:,3))).^2));

    %% --- 2) Standardized Opponent RGB ---
    Opp_F = reshape((reshape(rgbF,[],3) * M_opp.'), h, w, 3);
    Opp_H = reshape((reshape(rgbH,[],3) * M_opp.'), h, w, 3);

    % Compute differences using chromatic channels only (O1, O2)
    delta_all{2,k} = mat2gray(sqrt(sum((Opp_F(:,:,1:2)-Opp_H(:,:,1:2)).^2,3)));
    DoG_all{2,k}   = mat2gray(sqrt(sum((DoG(Opp_F(:,:,1:2))-DoG(Opp_H(:,:,1:2))).^2,3)));
end

%% --- Visualization ---
for s = 1:length(spaces)
    figure('Name',spaces{s},'Position',[50 50 1400 700]);
    for k = 1:3
        subplot(2,3,k); imagesc(delta_all{s,k}); axis image off; colorbar; colormap(jet);
        title([spaces{s} ' Δ – ' cases{k}{1}]);
        subplot(2,3,3+k); imagesc(DoG_all{s,k}); axis image off; colorbar; colormap(turbo);
        title([spaces{s} ' DoG – ' cases{k}{1}]);
    end
end












%% after vs before
%% --- Setup ---
sigma_small = 1;
sigma_large = 300;
DoG = @(X) imgaussfilt(X, sigma_small) - imgaussfilt(X, sigma_large);

RGB_before = painting_before.RGB_lin_img;
RGB_after  = painting_after.RGB_lin_img;

[h,w,~] = size(RGB_before);

%% --- Color spaces to compute differences ---
% 1) YCbCr (linear RGB version)
linear_rgb2ycbcr = @(rgb) cat(3, ...
    0.299*rgb(:,:,1) + 0.587*rgb(:,:,2) + 0.114*rgb(:,:,3), ...
    0.5*(rgb(:,:,3) - (0.299*rgb(:,:,1)+0.587*rgb(:,:,2)+0.114*rgb(:,:,3)))./(1-0.114+eps), ...
    0.5*(rgb(:,:,1) - (0.299*rgb(:,:,1)+0.587*rgb(:,:,2)+0.114*rgb(:,:,3)))./(1-0.299+eps) ...
);

% Preallocate
delta_all = cell(2,1);
DoG_all   = cell(2,1);

%% --- 1) YCbCr ---
YCbCr_before = linear_rgb2ycbcr(max(RGB_before,0));
YCbCr_after  = linear_rgb2ycbcr(max(RGB_after,0));

delta_all{1} = mat2gray(sqrt((YCbCr_before(:,:,2)-YCbCr_after(:,:,2)).^2 + ...
                              (YCbCr_before(:,:,3)-YCbCr_after(:,:,3)).^2));
DoG_all{1}   = mat2gray(sqrt((DoG(YCbCr_before(:,:,2))-DoG(YCbCr_after(:,:,2))).^2 + ...
                              (DoG(YCbCr_before(:,:,3))-DoG(YCbCr_after(:,:,3))).^2));

%% --- 2) Opponent RGB (standardized / orthonormal, De Valois & De Valois 1988) ---
M_opp = [ 1/sqrt(2)  -1/sqrt(2)   0;
           1/sqrt(6)   1/sqrt(6)  -2/sqrt(6);
           1/sqrt(3)   1/sqrt(3)   1/sqrt(3) ];

Opp_before = reshape(reshape(RGB_before,[],3) * M_opp.', h, w, 3);
Opp_after  = reshape(reshape(RGB_after, [],3) * M_opp.', h, w, 3);

% Compute Δ and DoG on the two chromatic channels (O1, O2)
delta_all{2} = mat2gray(sqrt((Opp_before(:,:,1)-Opp_after(:,:,1)).^2 + ...
                              (Opp_before(:,:,2)-Opp_after(:,:,2)).^2));

DoG_all{2}   = mat2gray(sqrt((DoG(Opp_before(:,:,1))-DoG(Opp_after(:,:,1))).^2 + ...
                              (DoG(Opp_before(:,:,2))-DoG(Opp_after(:,:,2))).^2));


%% --- Visualization ---
color_spaces = {'YCbCr','Opponent RGB'};

figure('Name','Δ and DoG maps before/after','Position',[100 100 900 700]);
for k = 1:2
    subplot(2,2,k); imagesc(delta_all{k}); axis image off; colorbar; colormap(turbo);
    title(['Δ - ' color_spaces{k}]);
    subplot(2,2,2+k); imagesc(DoG_all{k}); axis image off; colorbar; colormap(turbo);
    title(['DoG - ' color_spaces{k}]);
end







%%
%% --- Multi-scale DoG fusion for chromatic change detection ---
% Computes multi-scale DoG maps and fuses them (sum / max / L2)
% for Opponent RGB chroma (De Valois orthonormal) and YCbCr chroma.


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
[h, w, ~] = size(RGB_film);  % or whatever your base image variable is
min_dim = min(h, w);
% Adaptive large-scale sigma values based on image size
sigma_large_list = [min_dim/10, min_dim/5, min_dim/3];
nScales = numel(sigma_large_list);

% Use the same MSR cases setup you already have (film_case1/2/3 etc).
% Here I assume 'cases' variable exists as in your earlier code:
% cases = { {'No MSR', film_case1, hyper_case1}; ... };

% De Valois orthonormal opponent matrix (O1,O2,L)
M_opp = [ 1/sqrt(2), -1/sqrt(2), 0;
          1/sqrt(6),  1/sqrt(6), -2/sqrt(6);
          1/sqrt(3),  1/sqrt(3),  1/sqrt(3) ];

% helper DoG function
DoG_at = @(X, sL) imgaussfilt(X, sigma_small) - imgaussfilt(X, sL);

% preallocate
spaces = {'YCbCr','Opponent'};
delta_fused = cell(length(spaces),3); % sum, max, L2 fused map per case
DoG_fused = cell(length(spaces),3);



for k = 1:3
    rgbF = double(max(cases{k}{2}, 0));
    rgbH = double(max(cases{k}{3}, 0));
    [h,w,~] = size(rgbF);

    % --- Precompute Opponent and YCbCr chroma channels ---
    % YCbCr (linear-style chroma channels as before)
    linear_rgb2ycbcr = @(rgb) cat(3, ...
        0.299*rgb(:,:,1) + 0.587*rgb(:,:,2) + 0.114*rgb(:,:,3), ...
        0.5*(rgb(:,:,3) - (0.299*rgb(:,:,1)+0.587*rgb(:,:,2)+0.114*rgb(:,:,3)))./(1-0.114+eps), ...
        0.5*(rgb(:,:,1) - (0.299*rgb(:,:,1)+0.587*rgb(:,:,2)+0.114*rgb(:,:,3)))./(1-0.299+eps) ...
    );

    YCbCrF = linear_rgb2ycbcr(rgbF);
    YCbCrH = linear_rgb2ycbcr(rgbH);

    % Opponent (De Valois orthonormal) - reshape multiply for speed
    OppF = reshape((reshape(rgbF,[],3) * M_opp.'), h, w, 3);
    OppH = reshape((reshape(rgbH,[],3) * M_opp.'), h, w, 3);
    % Chromatic channels are the first two opponent dims (devalois ordering here)
    % (you can also use 2:3 depending on your chosen ordering; here O1,O2 are columns 1:2)

    % For each scale compute DoG differences on chroma channels and normalize
    DoG_maps_ycc = zeros(h,w,nScales);   % for YCbCr chroma magnitude per scale
    DoG_maps_opp = zeros(h,w,nScales);   % for Opponent chroma magnitude per scale

    for si = 1:nScales
        sL = sigma_large_list(si);

        % YCbCr chroma channels (indices 2 and 3)
        dCb = DoG_at(YCbCrF(:,:,2), sL) - DoG_at(YCbCrH(:,:,2), sL);
        dCr = DoG_at(YCbCrF(:,:,3), sL) - DoG_at(YCbCrH(:,:,3), sL);
        mag_ycc = sqrt(dCb.^2 + dCr.^2);

        % Opponent chroma channels (use O1 and O2). depending on ordering:
        % Using columns 1 and 2 as chroma (O1, O2)
        dO1 = DoG_at(OppF(:,:,1), sL) - DoG_at(OppH(:,:,1), sL);
        dO2 = DoG_at(OppF(:,:,2), sL) - DoG_at(OppH(:,:,2), sL);
        mag_opp = sqrt(dO1.^2 + dO2.^2);

        DoG_maps_ycc(:,:,si) = mag_ycc;
        DoG_maps_opp(:,:,si) = mag_opp;
    end

    % --- Robust normalization per scale to suppress noise dominance ---
    % Options: use MAD or high-percentile. I'll use MAD here.
    normalize_by_mad = @(M) (M ./ (repmat(mad(reshape(M,[],1),1), size(M)) + eps));
    % But mad returns a scalar for the whole map; we want per-scale scalar
    % for si = 1:nScales
    %     v = DoG_maps_ycc(:,:,si); DoG_maps_ycc(:,:,si) = v ./ (mad(v(:),1) + eps);
    %     v = DoG_maps_opp(:,:,si); DoG_maps_opp(:,:,si) = v ./ (mad(v(:),1) + eps);
    % end

    % --- Fusion strategies: sum, max, L2 (rss) ---
    fused_sum_ycc = sum(DoG_maps_ycc, 3);
    fused_max_ycc = max(DoG_maps_ycc, [], 3);
    fused_l2_ycc  = sqrt(sum(DoG_maps_ycc.^2, 3));

    fused_sum_opp = sum(DoG_maps_opp, 3);
    fused_max_opp = max(DoG_maps_opp, [], 3);
    fused_l2_opp  = sqrt(sum(DoG_maps_opp.^2, 3));

    % Normalize for display with mat2gray
    delta_fused{1,k} = mat2gray(fused_sum_ycc); % YCbCr sum
    DoG_fused{1,k}   = mat2gray(fused_l2_ycc);  % YCbCr L2 (choose one for DoG visualization)

    delta_fused{2,k} = mat2gray(fused_sum_opp); % Opponent sum
    DoG_fused{2,k}   = mat2gray(fused_l2_opp);

    % also keep max variants if you want to inspect
    % delta_fused_max_ycc{k} = mat2gray(fused_max_ycc); ...
end

%% --- Visualization (sum-fused and L2-DoG fused) ---
figure('Name','Multi-scale fused chroma (YCbCr sum / Opp sum)','Position',[100 100 1200 600]);
for k = 1:3
    subplot(2,3,k); imagesc(delta_fused{1,k}); axis image off; colorbar; colormap(turbo);
    title(['YCbCr fused (sum) – ' cases{k}{1}]);
    subplot(2,3,3+k); imagesc(DoG_fused{1,k}); axis image off; colorbar; colormap(turbo);
    title(['YCbCr fused (L2 DoG) – ' cases{k}{1}]);
end

figure('Name','Multi-scale fused chroma (Opponent sum / Opp L2)','Position',[100 100 1200 600]);
for k = 1:3
    subplot(2,3,k); imagesc(delta_fused{2,k}); axis image off; colorbar; colormap(turbo);
    title(['Opponent fused (sum) – ' cases{k}{1}]);
    subplot(2,3,3+k); imagesc(DoG_fused{2,k}); axis image off; colorbar; colormap(turbo);
    title(['Opponent fused (L2 DoG) – ' cases{k}{1}]);
end





%%









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







%% XYZ
%% --- Parameters ---
scales = [15, 80, 250];  % MSR scales
sigma_small = 1;
sigma_large = 300;
DoG = @(X) imgaussfilt(X, sigma_small) - imgaussfilt(X, sigma_large);

%% --- Load XYZ images ---
XYZ_film  = film_data.XYZ_img;       % HxWx3
XYZ_hyper = painting_after.XYZ_img;  % HxWx3
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




