clear; close all;

% --- Paths
path_before = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat';
path_after  = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat';
path_film   = '/home/oem/eliza/data/xyz_lab_rgb/film/yoda_halogen_fuji_exp0.mat';
% path_before = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat';
% path_after  = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat';
% path_film   = '/home/oem/eliza/data/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat';

outputDir = '/home/oem/eliza/masters-thesis/results/plots/dog';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
[~, filmName, ~] = fileparts(path_film);

% --- Load
painting_before = load(path_before);
painting_after  = load(path_after);
film_data       = load(path_film);

% --- Lab images
Lab_before = painting_before.Lab_img;
Lab_after  = painting_after.Lab_img;
Lab_film   = film_data.Lab_img;

% --- Multi-scale DoG on Lab channels (suppress global differences)
sigma_pairs = [1 10; 1.5 85; 2 780];
Lab_film_dog  = multiscale_dog_filter(Lab_film, sigma_pairs);
Lab_after_dog = multiscale_dog_filter(Lab_after, sigma_pairs);

% --- ΔE2000 comparison (DoG-processed)
dE_map_dog    = compute_deltaE2000(Lab_film_dog, Lab_after);
dE_map_direct = compute_deltaE2000(Lab_before, Lab_after);  % ground truth unmodified

% --- Thresholded masks
fixed_threshold = 6;
mask_dog    = dE_map_dog    >= fixed_threshold;
mask_direct = dE_map_direct >= fixed_threshold;

% --- Display sRGB
film_srgb  = xyz2rgb(film_data.XYZ_img ./ 100, 'WhitePoint', 'd50');
after_srgb = xyz2rgb(painting_after.XYZ_img ./ 100, 'WhitePoint', 'd50');





%% --- Figure  Change Detection Masks ---
h2 = figure('Name', 'Change Detection Masks', 'Units', 'normalized', 'Position', [0.05 0.1 0.75 0.6]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile; imshow(mask_dog);    title({'Multi-Scale DoG ΔE00','Film vs After'});
nexttile; imshow(mask_direct); title({'Ground Truth','Before vs After ΔE ≥ 6'});

% export figure 2
file2 = fullfile(outputDir, sprintf('%s_Change_Detection_Masks.png', filmName));
% exportgraphics(h2, file2, 'Resolution', 300);




% --- Function: Multi-scale DoG
function Lab_dog = multiscale_dog_filter(Lab_img, sigma_pairs)
    Lab_dog = zeros(size(Lab_img));
    for c = 1:3
        acc = zeros(size(Lab_img(:,:,1)));
        for i = 1:size(sigma_pairs,1)
            blur1 = imgaussfilt(Lab_img(:,:,c), sigma_pairs(i,1));
            blur2 = imgaussfilt(Lab_img(:,:,c), sigma_pairs(i,2));
            acc = acc + (blur1 - blur2);
        end
        acc = acc / size(sigma_pairs,1);  % Average over scales
        Lab_dog(:,:,c) = acc;
    end
end

% --- Function: Compute ΔE2000
function dE_map = compute_deltaE2000(Lab1, Lab2)
    sz = size(Lab1);
    Lab1 = reshape(Lab1, [], 3);
    Lab2 = reshape(Lab2, [], 3);
    dE = deltaE2000(Lab1, Lab2);
    dE_map = reshape(dE, sz(1), sz(2));
end
