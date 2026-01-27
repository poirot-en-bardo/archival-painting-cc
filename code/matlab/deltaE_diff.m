%% User-defined image paths
film_path = '/Volumes/School/Thesis/thesis-repo/results/plots/cc/target_based/yoda_halogen_fuji_exp0_/yoda_halogen_fuji_exp0_corrected_rgb_srgb.tif';
hsi_before_path = '/Volumes/School/Thesis/thesis-repo/results/plots/srgb/yoda_reflectance_before_xyz_srgb.png';
hsi_after_path = '/Volumes/School/Thesis/thesis-repo/results/plots/srgb/yoda_reflectance_after_reg_xyz_srgb.png';
[~, film_base_name, ~] = fileparts(film_path);
base_output_dir = '/Volumes/School/Thesis/thesis-repo/results/plots/cc/target_based/residual';
outputDir = fullfile(base_output_dir);
%% Load images
film_img = im2double(imread(film_path));
hsi_before = im2double(imread(hsi_before_path));
hsi_after = im2double(imread(hsi_after_path));

%% Convert sRGB to CIELAB
cform = makecform('srgb2lab');
film_lab = applycform(film_img, cform);
hsi_before_lab = applycform(hsi_before, cform);
hsi_after_lab = applycform(hsi_after, cform);

%% Reshape to N x 3 for deltaE function
[nrows, ncols, ~] = size(film_lab);
film_lab_reshaped = reshape(film_lab, [], 3);
hsi_before_reshaped = reshape(hsi_before_lab, [], 3);
hsi_after_reshaped = reshape(hsi_after_lab, [], 3);

%% Compute deltaE maps
DE_film_hsi = deltaE2000(film_lab_reshaped, hsi_after_reshaped);
DE_hsi_before_after = deltaE2000(hsi_before_reshaped, hsi_after_reshaped);

%% Reshape deltaE back to image
DE_film_hsi_img = reshape(DE_film_hsi, nrows, ncols);
DE_hsi_before_after_img = reshape(DE_hsi_before_after, nrows, ncols);

%% Compute difference map
DE_diff = DE_film_hsi_img - DE_hsi_before_after_img;

%% Display difference map
% Parameters
figSize = [0.1, 0.2, 1, 0.8];  % [left, bottom, width, height]
clim_range = [0 20];
fontSize = 20;
labelSize = 30;          % ΔE scale

% Create figure
fig2 = figure('Units', 'normalized', 'Position', figSize);

% Display the ΔE difference
imagesc(DE_diff);
axis image off;  % keep axes off
colormap parula;
clim(clim_range);

% Title
title('Correction Error', ...
      'FontSize', 15, 'FontWeight', 'bold');

% Colorbar
c = colorbar;
c.Label.String = 'ΔE_{00}';
c.Label.FontSize = labelSize;
c.Label.FontWeight = 'bold';

% Set axes font
set(gca, 'FontSize', labelSize, 'FontWeight', 'bold');

% Save figure
exportgraphics(fig2, fullfile(outputDir, [film_base_name '_deltaE_residual.png']), ...
               'Resolution', 300, 'BackgroundColor', 'none', 'ContentType', 'image');


%% Save the figure
saveas(gcf, 'deltaE_difference.png');
