%% --- MTT-based change detection with MSR cases ---
clear; close all;

%% --- Paths ---
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

%% --- Load data ---
before = load(path_before);
after  = load(path_after);
film   = load(path_film);

% Linear RGB in [0,1]
RGB_before = before.RGB_lin_img;
RGB_after  = after.RGB_lin_img;
RGB_film   = film.RGB_lin_img;

[h,w,~] = size(RGB_film);
min_dim = min(h,w);

%% --- Parameters ---
eps_val = 1e-9;
DoG_small = 1;
sigma_large_list = [min_dim/10, min_dim/5, min_dim/3];
normalizeLuma = true;
clipChromRange = true;
msr_scales = [15 80 250]; % MSR scales

DoG_at = @(X,sL) imgaussfilt(X, DoG_small) - imgaussfilt(X, sL);

%% --- Helper Functions ---

% Piecewise angular warp (oRGB paper)
function theta_o = piecewise_theta_o(theta)
    theta2 = mod(theta + 2*pi, 2*pi);
    theta_o2 = zeros(size(theta2));
    half_mask = theta2 <= pi;
    t_half = theta2(half_mask);
    theta_o2(half_mask) = (3/2) * t_half;
    idx_p2 = half_mask & (theta2 > pi/3);
    t2 = theta2(idx_p2);
    theta_o2(idx_p2) = pi/2 + (3/4)*(t2 - pi/3);
    idx_reflect = theta2 > pi;
    t_ref = 2*pi - theta2(idx_reflect);
    theta_ref_out = zeros(size(t_ref));
    mask_r1 = t_ref <= (pi/3);
    theta_ref_out(mask_r1) = (3/2) * t_ref(mask_r1);
    theta_ref_out(~mask_r1) = pi/2 + (3/4)*(t_ref(~mask_r1)-pi/3);
    theta_o2(idx_reflect) = -theta_ref_out;
    theta_o = theta_o2;
    theta_o(theta_o>pi) = theta_o(theta_o>pi) - 2*pi;
end

% Compute oRGB chroma channels
function [C1,C2] = compute_oRGB_chroma(R,G,B,normalizeLuma,clipChromRange)
    eps_val = 1e-9;
    L = 0.299*R + 0.587*G + 0.114*B;
    yb = 0.5*(R+G)-B;
    rg = 0.866*(R-G);
    if normalizeLuma
        Lmin = min(L(:)); Lmax = max(L(:));
        L = (L-Lmin)/(Lmax-Lmin+eps_val);
    end
    theta = atan2(rg,yb);
    theta_o = piecewise_theta_o(theta);
    delta = theta_o - theta;
    C1 = cos(delta).*yb - sin(delta).*rg;
    C2 = sin(delta).*yb + cos(delta).*rg;
    if clipChromRange
        maxMag = max([max(abs(C1(:))), max(abs(C2(:)))]);
        if maxMag>0
            C1 = C1 / (maxMag+eps_val);
            C2 = C2 / (maxMag+eps_val);
        end
    end
end

% Multi-Scale Retinex (MSR) function
function img_msr = MSR_nl(img, scales)
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
        eps_val = 1e-8;
        msr_channel = (msr_channel - min(msr_channel(:))) / (max(msr_channel(:)) - min(msr_channel(:)) + eps_val);
        img_msr(:,:,c) = msr_channel;
    end
end

% Compute Δ-chroma and DoG fused maps
function [delta_chroma,DoG_fused] = change_maps(C1A,C2A,C1B,C2B,sigma_large_list,DoG_at)
    delta_chroma = sqrt((C1A-C1B).^2 + (C2A-C2B).^2);
    delta_chroma = mat2gray(delta_chroma);
    nScales = numel(sigma_large_list);
    [h,w] = size(C1A);
    DoG_maps = zeros(h,w,nScales);
    for si=1:nScales
        sL = sigma_large_list(si);
        d1 = DoG_at(C1A,sL) - DoG_at(C1B,sL);
        d2 = DoG_at(C2A,sL) - DoG_at(C2B,sL);
        DoG_maps(:,:,si) = sqrt(d1.^2 + d2.^2);
    end
    DoG_fused = mat2gray(mean(DoG_maps,3));
end

%% --- Build MSR cases (only Film vs After has MSR) ---

film_MSR =  MSR_nl(RGB_film, msr_scales);
after_MSR =  MSR_nl(RGB_film, msr_scales);
film_MSR_gamma = gamma_encode(film_MSR);
after_MSR_gamma = gamma_encode(after_MSR);
cases = { ...
    {'No MSR', RGB_film, RGB_after, RGB_before}, ...
    {'Film MSR', film_MSR_gamma, RGB_after, RGB_before}, ...
    {'Both MSR', film_MSR_gamma, after_MSR_gamma, RGB_before} ...
};

%% --- Store results for plotting ---
delta_chroma_FA_cases = cell(1,numel(cases));
DoG_FA_cases = cell(1,numel(cases));

% Only need Before vs After once
delta_chroma_BA = [];
DoG_BA = [];

%% --- Process each case ---
for k = 1:numel(cases)
    case_name = cases{k}{1};
    imgF = cases{k}{2};  % Film (or MSR)
    imgA = cases{k}{3};  % After (MSR only if Both MSR)
    imgB = cases{k}{4};  % Before always original linear RGB
    
    % Chroma channels
    [C1_F,C2_F] = compute_oRGB_chroma(imgF(:,:,1),imgF(:,:,2),imgF(:,:,3),normalizeLuma,clipChromRange);
    [C1_A,C2_A] = compute_oRGB_chroma(imgA(:,:,1),imgA(:,:,2),imgA(:,:,3),normalizeLuma,clipChromRange);
    [C1_B,C2_B] = compute_oRGB_chroma(imgB(:,:,1),imgB(:,:,2),imgB(:,:,3),normalizeLuma,clipChromRange);
    
    % Change maps
    [delta_chroma_FA,DoG_FA] = change_maps(C1_F,C2_F,C1_A,C2_A,sigma_large_list,DoG_at);
    [delta_chroma_BA_tmp,DoG_BA_tmp] = change_maps(C1_B,C2_B,C1_A,C2_A,sigma_large_list,DoG_at);
    
    % Store Film vs After results
    delta_chroma_FA_cases{k} = delta_chroma_FA;
    DoG_FA_cases{k} = DoG_FA;
    
    % Store Before vs After only once (they are the same for all MSR cases)
    if isempty(delta_chroma_BA)
        delta_chroma_BA = delta_chroma_BA_tmp;
        DoG_BA = DoG_BA_tmp;
    end
end

%% --- Plot Before vs After ---
figure('Name','Before vs After','Position',[100 100 1000 500]);
subplot(1,2,1);
imagesc(delta_chroma_BA); axis image off; colorbar; colormap(turbo); title('Δ Chroma: Before vs After');

subplot(1,2,2);
imagesc(DoG_BA); axis image off; colorbar; colormap(turbo); title('DoG Fused: Before vs After');

%% --- Plot Film vs After MSR cases ---
figure('Name','Film vs After MSR Cases','Position',[100 100 1400 800]);
msr_case_names = {'No MSR','Film MSR','Both MSR'};

for k = 1:3
    % Δ Chroma
    subplot(2,3,k);
    imagesc(delta_chroma_FA_cases{k}); axis image off; colorbar; colormap(turbo);
    title(['Δ Chroma: ' msr_case_names{k}]);
    
    % DoG
    subplot(2,3,3+k);
    imagesc(DoG_FA_cases{k}); axis image off; colorbar; colormap(turbo);
    title(['DoG: ' msr_case_names{k}]);
end


function img_gamma = gamma_encode(img_lin, gamma)
%GAMMA_ENCODE Apply gamma encoding to a linear RGB image
%   img_lin: linear RGB image in [0,1]
%   gamma: gamma value (default 2.2 if not specified)
%
%   Returns: gamma-encoded image in [0,1]

if nargin < 2
    gamma = 2.2;
end

% Ensure values are in [0,1]
img_lin = max(min(img_lin, 1), 0);

% Apply gamma encoding
img_gamma = img_lin .^ (1/gamma);

end
