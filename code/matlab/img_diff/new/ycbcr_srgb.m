
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


%%




%% --- Setup ---
scales = [15, 80, 250]; % MSR scales
sigma_small = 1;
sigma_large = 300;
DoG = @(X) imgaussfilt(X, sigma_small) - imgaussfilt(X, sigma_large);

RGB_film  = film_data.RGB_lin_img;
RGB_hyper = painting_after.RGB_lin_img;

[h,w,~] = size(RGB_film);

%% --- Prepare 3 MSR cases ---
XYZ_film = film_data.XYZ_img;
XYZ_hyper = painting_after.XYZ_img;

% RGB_film = xyz2rgb(XYZ_film,"ColorSpace","linear-rgb", "OutputType","double","WhitePoint","d65");
% RGB_hyper = xyz2rgb(XYZ_hyper,"ColorSpace","linear-rgb", "OutputType","double", "WhitePoint","d65");
RGB_film  = film_data.RGB_lin_img;
RGB_hyper = painting_after.RGB_lin_img;

film_case1 = RGB_film;
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


sanitize_rgb = @(R) sanitize_rgb_fn(R);

function R = sanitize_rgb_fn(R_in)
    % convert to double
    R = double(R_in);
    % remove tiny imaginary parts if any
    if ~isreal(R)
        R = real(R);
    end
    % replace NaN/Inf
    R(~isfinite(R)) = 0;
    % clamp to [0,1] (rgb2ycbcr expects doubles in [0,1] for double input)
    R = min(max(R, 0), 1);
end

%% --- Initialize storage ---
spaces = {'YCbCr'};
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

    % rgbF = sanitize_rgb(cases{k}{2});
    % rgbH = sanitize_rgb(cases{k}{3});

    %% --- 1) YCbCr ---
    YCbCr_F = rgb2ycbcr(rgbF);
    YCbCr_H = rgb2ycbcr(rgbH);
    delta_all{1,k} = mat2gray(sqrt((YCbCr_F(:,:,2)-YCbCr_H(:,:,2)).^2 + ...
                                    (YCbCr_F(:,:,3)-YCbCr_H(:,:,3)).^2));
    DoG_all{1,k}   = mat2gray(sqrt((DoG(YCbCr_F(:,:,2))-DoG(YCbCr_H(:,:,2))).^2 + ...
                                    (DoG(YCbCr_F(:,:,3))-DoG(YCbCr_H(:,:,3))).^2));

   
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






%% ----------------------------------------------------------
%  LOAD FILES
%% ----------------------------------------------------------



XYZ_before = painting_before.XYZ_img;
XYZ_after  = painting_after.XYZ_img;

%% ----------------------------------------------------------
%  Convert XYZ → linear RGB (same as your pipeline)
%% ----------------------------------------------------------

RGB_before = xyz2rgb(XYZ_before, ...
    "ColorSpace","linear-rgb", ...
    "OutputType","double", ...
    "WhitePoint","d65");

RGB_after = xyz2rgb(XYZ_after, ...
    "ColorSpace","linear-rgb", ...
    "OutputType","double", ...
    "WhitePoint","d65");

%% Clamp to [0,1] because rgb2ycbcr requires double in [0,1]
RGB_before = max(RGB_before,0);
RGB_after  = max(RGB_after,0);

RGB_before = painting_before.RGB_lin_img;
RGB_after = painting_after.RGB_lin_img;

[h,w,~] = size(RGB_before);

%% ----------------------------------------------------------
%  DoG SETUP (same as yours)
%% ----------------------------------------------------------

sigma_small = 1;
sigma_large = 300;

DoG = @(X) imgaussfilt(X,sigma_small) - imgaussfilt(X,sigma_large);

%% ----------------------------------------------------------
%  Convert to YCbCr (MATLAB built-in, same as your code)
%% ----------------------------------------------------------

YCbCr_B = rgb2ycbcr(RGB_before);
YCbCr_A = rgb2ycbcr(RGB_after);

%% ----------------------------------------------------------
%  DELTA (same structure as your code)
%% ----------------------------------------------------------

delta_YCbCr = mat2gray( sqrt( ...
    (YCbCr_B(:,:,2) - YCbCr_A(:,:,2)).^2 + ...
    (YCbCr_B(:,:,3) - YCbCr_A(:,:,3)).^2 ) );

%% ----------------------------------------------------------
%  DoG Δ (same structure as your code)
%% ----------------------------------------------------------

DoG_YCbCr = mat2gray( sqrt( ...
    (DoG(YCbCr_B(:,:,2)) - DoG(YCbCr_A(:,:,2))).^2 + ...
    (DoG(YCbCr_B(:,:,3)) - DoG(YCbCr_A(:,:,3))).^2 ) );

%% ----------------------------------------------------------
%  Visualization (identical style to your script)
%% ----------------------------------------------------------

figure('Name','YCbCr Before vs After','Position',[50 50 1200 600]);

subplot(1,2,1);
imagesc(delta_YCbCr); axis image off; colorbar; colormap(jet);
title('Δ (YCbCr chromatic)');

subplot(1,2,2);
imagesc(DoG_YCbCr); axis image off; colorbar; colormap(turbo);
title('DoG Δ (YCbCr chromatic)');
