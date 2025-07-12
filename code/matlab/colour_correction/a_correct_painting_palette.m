clear; close all;

%% === File Paths ===
painting_after_path = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat';
painting_before_path = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat';
% painting_after_path = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat';
% painting_before_path = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat';
film_path = '/home/oem/eliza/data/xyz_lab_rgb/film/yoda_led_kodak_exp0.mat';
cluster_rgb_file = '/home/oem/eliza/data/xyz_lab_rgb/clustered/filtered/yoda_led_kodak_exp0/yoda_led_kodak_exp0_cluster.mat';
painting_after_palette_file = '/home/oem/eliza/data/xyz_lab_rgb/clustered/yoda_painting_after_cluster_rgb.mat';



%% === Load Clustered Data ===
film_cluster_rgb = load(cluster_rgb_file).rgb_film_clustered;
painting_after_cluster_rgb = load(painting_after_palette_file).rgb_after_clustered;

painting_after = load(painting_after_path);
painting_before = load(painting_before_path);
film_data = load(film_path);

%%

% Extract base name from film image path
[~, film_base_name, ~] = fileparts(film_path);

% Define output directory using the film base name
base_output_dir = '/home/oem/eliza/masters-thesis/results/plots/cc/palette';
outputDir = fullfile(base_output_dir, film_base_name);

% Create the directory if it doesn't exist
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
%%
cc_rgb_ref = painting_after_cluster_rgb;
cc_rgb_film = film_cluster_rgb;
painting_after_rgb = painting_after.RGB_img;
painting_before_rgb = painting_before.RGB_img;
film_rgb = film_data.RGB_img;

%% Colour correction
% Prepare training data from colour checker
X_train = root_poly_features(cc_rgb_film);   % Input: film RGB patches
Y_train = cc_rgb_ref;                        % Target: spectral RGB patches

lambda = 0.001;
coeffs = pinv(X_train' * X_train + lambda * eye(size(X_train,2))) * X_train' * Y_train;

%% Apply correction to full film image
[h, w, ~] = size(film_rgb);
X_film = root_poly_features(reshape(film_rgb, [], 3));  % Expand film image RGB
corrected_rgb = X_film * coeffs;                        % Apply regression
corrected_rgb = min(max(corrected_rgb, 0), 1);          % Clamp to [0, 1]
corrected_rgb_img = reshape(corrected_rgb, h, w, 3);

%% Convert both corrected and reference to Lab via XYZ
% ===== Convert corrected and reference images to XYZ =====
xyz_corr       = prophoto2xyz(corrected_rgb, true);  % Film after correction
xyz_ref_after  = prophoto2xyz(reshape(painting_after_rgb, [], 3), true);  % Spectral (after ageing)
xyz_ref_before = prophoto2xyz(reshape(painting_before_rgb, [], 3), true); % Spectral (before ageing)

% ===== Convert XYZ to Lab =====
lab_corr       = xyz2lab_custom(xyz_corr);
lab_ref_after  = xyz2lab_custom(xyz_ref_after);
lab_ref_before = xyz2lab_custom(xyz_ref_before);

% ===== Compute ΔE2000 =====
deltaE_after  = deltaE2000(lab_ref_after, lab_corr);
deltaE_before = deltaE2000(lab_ref_before, lab_corr);

% ===== Reshape to image form =====
deltaE_after_img  = reshape(deltaE_after, h, w);
deltaE_before_img = reshape(deltaE_before, h, w);

% Apply the masks to ΔE maps
mask_after  = painting_after.spec_mask;   % true = specular pixel
mask_before = painting_before.spec_mask;

deltaE_after_img(mask_after)   = 0;
deltaE_before_img(mask_before) = 0;

deltaE_after_thresholded = deltaE_after_img;
deltaE_after_thresholded(deltaE_after_thresholded < 10) = 0;

% --- Compute ΔE between painting before and after ageing ---
lab_before = reshape(painting_before.Lab_img, [], 3);
lab_after  = reshape(painting_after.Lab_img, [], 3);

deltaE_painting = deltaE2000(lab_before, lab_after);
deltaE_painting_img = reshape(deltaE_painting, h, w);

% --- Apply both masks ---
combined_mask = painting_before.spec_mask | painting_after.spec_mask;
deltaE_painting_img(combined_mask) = 0;

%% Visualisation
% Convert all RGB images from ProPhoto to sRGB for display
film_rgb_srgb          = xyz2rgb(film_data.XYZ_img ./100, 'ColorSpace', 'srgb', 'WhitePoint','d50');
corrected_rgb_srgb     = xyz2rgb(prophoto2xyz(corrected_rgb, true) ./100, 'ColorSpace', 'srgb', 'WhitePoint','d50');
painting_before_srgb   = xyz2rgb(painting_before.XYZ_img ./100, 'ColorSpace', 'srgb', 'WhitePoint', 'd50');
painting_after_srgb    = xyz2rgb(painting_after.XYZ_img ./100, 'ColorSpace', 'srgb', 'WhitePoint','d50');

% Reshape to image form
film_rgb_srgb_img        = reshape(film_rgb_srgb, h, w, 3);
corrected_rgb_srgb_img   = reshape(corrected_rgb_srgb, h, w, 3);
painting_before_srgb_img = reshape(painting_before_srgb, h, w, 3);
painting_after_srgb_img  = reshape(painting_after_srgb, h, w, 3);
%%
% Visualise all in sRGB
fig_vis = figure('Units', 'normalized', 'Position', [0.05, 0.05, 1.2, 1]);
tiledlayout(fig_vis, 1, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

fontSize = 20;
fontWeight = 'bold';

nexttile;
imshow(film_rgb_srgb_img);
title('Film Before Correction', 'FontSize', fontSize, 'FontWeight', fontWeight);

nexttile;
imshow(corrected_rgb_srgb_img);
title('Film After Correction', 'FontSize', fontSize, 'FontWeight', fontWeight);

nexttile;
imshow(painting_before_srgb_img);
title('Painting (Before Ageing)', 'FontSize', fontSize, 'FontWeight', fontWeight);
% 
% nexttile;
% imshow(painting_after_srgb_img);
% title('Painting (After Ageing)', 'FontSize', fontSize, 'FontWeight', fontWeight);

%%
% Save figures
fig_vis_path = fullfile(outputDir, [film_base_name '_visualisation.png']);
% saveas(fig_vis, fig_vis_path);
exportgraphics(fig_vis, fig_vis_path, 'Resolution', 300, 'BackgroundColor', 'none', 'ContentType', 'image');

fprintf('Saved visualisation figure to: %s\n', fig_vis_path);



%% Delta E maps
% Define common figure position and size
% figSize = [0.1, 0.2, 1, 0.8];  % [left, bottom, width, height]
% clim_range = [0 20];
% fontSize = 20;
% labelSize = 20;
% 
% % --- 1. ΔE map: corrected vs. HSI before ageing ---
% fig1 = figure('Units', 'normalized', 'Position', figSize);
% imagesc(deltaE_before_img);
% axis image off;
% title('Corrected film vs. HSI before ageing', 'FontSize', fontSize, 'FontWeight', 'bold');
% c = colorbar;
% clim(clim_range);
% c.Label.String = 'ΔE_{00}';
% c.Label.FontSize = labelSize;
% c.Label.FontWeight = 'bold';
% set(gca, 'FontSize', labelSize, 'FontWeight', 'bold');
% % saveas(fig1, fullfile(outputDir, [film_base_name '_deltaE_before.png']));
% exportgraphics(fig1, fullfile(outputDir, [film_base_name '_deltaE_before.png']), 'Resolution', 300, 'BackgroundColor', 'none', 'ContentType', 'image');
% 
% 
% % --- 2. ΔE map: corrected vs. HSI after ageing ---
% fig2 = figure('Units', 'normalized', 'Position', figSize);
% imagesc(deltaE_after_img);
% axis image off;
% title('Corrected film vs. HSI after ageing', 'FontSize', fontSize, 'FontWeight', 'bold');
% c = colorbar;
% clim(clim_range);
% c.Label.String = 'ΔE_{00}';
% c.Label.FontSize = labelSize;
% c.Label.FontWeight = 'bold';
% set(gca, 'FontSize', labelSize, 'FontWeight', 'bold');
% % saveas(fig2, fullfile(outputDir, [film_base_name '_deltaE_after.png']));
% exportgraphics(fig2, fullfile(outputDir, [film_base_name '_deltaE_after.png']), 'Resolution', 300, 'BackgroundColor', 'none', 'ContentType', 'image');
% 
% 
% % --- 3. ΔE map: HSI before vs. after ageing ---
% fig3 = figure('Units', 'normalized', 'Position', figSize);
% imagesc(deltaE_painting_img);
% axis image off;
% title('HSI before vs. after ageing', 'FontSize', fontSize, 'FontWeight', 'bold');
% c = colorbar;
% clim(clim_range);
% c.Label.String = 'ΔE_{00}';
% c.Label.FontSize = labelSize;
% c.Label.FontWeight = 'bold';
% set(gca, 'FontSize', labelSize, 'FontWeight', 'bold');
% % saveas(fig3, fullfile(outputDir, [film_base_name '_deltaE_painting_before_after.png']));
% exportgraphics(fig3, fullfile(outputDir, [film_base_name '_deltaE_painting_before_after.png']), 'Resolution', 300, 'BackgroundColor', 'none', 'ContentType', 'image');
% 
% 
% % --- 4. ΔE map: thresholded ΔE > 10 ---
% fig4 = figure('Units', 'normalized', 'Position', figSize);
% imagesc(deltaE_after_thresholded);
% axis image off;
% title('Corrected film vs. after ageing (ΔE > 10)', 'FontSize', fontSize, 'FontWeight', 'bold');
% c = colorbar;
% clim([10 20]);
% c.Label.String = 'ΔE_{00}';
% c.Label.FontSize = labelSize;
% c.Label.FontWeight = 'bold';
% set(gca, 'FontSize', labelSize, 'FontWeight', 'bold');
% % saveas(fig4, fullfile(outputDir, [film_base_name '_deltaE_after_thresh10.png']));
% exportgraphics(fig4, fullfile(outputDir, [film_base_name '_deltaE_after_thresh10.png']), 'Resolution', 300, 'BackgroundColor', 'none', 'ContentType', 'image');


%% === ΔE Maps (Combined 2x2 Layout) ===
fig_combined = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 1.05]);
tiledlayout(fig_combined, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

clim_range = [0 20];
fontSize = 18;
labelSize = 18;
titleYOffset = 1.03;  % move title up slightly

% 1. Corrected vs Before
nexttile;
imagesc(deltaE_before_img); axis image off;
t = title('Corrected film vs. HSI before ageing', 'FontSize', fontSize, 'FontWeight', 'bold');
t.Units = 'normalized'; t.Position(2) = titleYOffset;
c = colorbar; clim(clim_range);
c.Label.String = '\DeltaE_{00}'; c.Label.FontSize = labelSize; c.Label.FontWeight = 'bold';
set(c, 'FontSize', labelSize, 'FontWeight', 'bold');

% 2. Corrected vs After
nexttile;
imagesc(deltaE_after_img); axis image off;
t = title('Corrected film vs. HSI after ageing', 'FontSize', fontSize, 'FontWeight', 'bold');
t.Units = 'normalized'; t.Position(2) = titleYOffset;
c = colorbar; clim(clim_range);
c.Label.String = '\DeltaE_{00}'; c.Label.FontSize = labelSize; c.Label.FontWeight = 'bold';
set(c, 'FontSize', labelSize, 'FontWeight', 'bold');

% 3. HSI Before vs After
nexttile;
imagesc(deltaE_painting_img); axis image off;
t = title('HSI before vs. HSI after ageing (GT)', 'FontSize', fontSize, 'FontWeight', 'bold');
t.Units = 'normalized'; t.Position(2) = titleYOffset;
c = colorbar; clim(clim_range);
c.Label.String = '\DeltaE_{00}'; c.Label.FontSize = labelSize; c.Label.FontWeight = 'bold';
set(c, 'FontSize', labelSize, 'FontWeight', 'bold');

% 4. Thresholded ΔE
nexttile;
imagesc(deltaE_after_thresholded); axis image off;
t = title('Corrected film vs. HSI after ageing (\DeltaE > 10)', 'FontSize', fontSize, 'FontWeight', 'bold');
t.Units = 'normalized'; t.Position(2) = titleYOffset;
c = colorbar; clim([10 20]);
c.Label.String = '\DeltaE_{00}'; c.Label.FontSize = labelSize; c.Label.FontWeight = 'bold';
set(c, 'FontSize', labelSize, 'FontWeight', 'bold');

drawnow; pause(0.5);
exportgraphics(fig_combined, fullfile(outputDir, [film_base_name '_deltaE_maps_combined.png']), 'Resolution', 300);

%% Delta E histograms
% %% Plot histograms of ΔE maps
% fig_hist = figure('Units', 'normalized', 'Position', [0.1, 0.2, 1, 0.6]);
% tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
% 
% edges = 0:1:40;  % Bin edges for ΔE values
% labelFontSize = 14;
% tickFontSize = 15;
% 
% % 1. Histogram: ΔE before
% nexttile;
% histogram(deltaE_before_img(:), edges, 'FaceColor', [0.2 0.6 0.8]);
% title('Corrected Film vs. Before Ageing', 'FontSize', 18, 'FontWeight', 'bold');
% xlabel('ΔE_{00}', 'FontSize', labelFontSize, 'FontWeight', 'bold');
% ylabel('Pixel Count', 'FontSize', labelFontSize, 'FontWeight', 'bold');
% xlim([0 40]); grid on;
% set(gca, 'FontSize', tickFontSize, 'FontWeight', 'bold');
% 
% % 2. Histogram: ΔE after
% nexttile;
% histogram(deltaE_after_img(:), edges, 'FaceColor', [0.2 0.8 0.4]);
% title('Corrected Film vs. After Ageing', 'FontSize', 18, 'FontWeight', 'bold');
% xlabel('ΔE_{00}', 'FontSize', labelFontSize, 'FontWeight', 'bold');
% ylabel('Pixel Count', 'FontSize', labelFontSize, 'FontWeight', 'bold');
% xlim([0 40]); grid on;
% set(gca, 'FontSize', tickFontSize, 'FontWeight', 'bold');
% 
% % 3. Histogram: Painting Before vs. After
% nexttile;
% histogram(deltaE_painting_img(:), edges, 'FaceColor', [0.9 0.5 0.2]);
% title('HSI Painting Before vs. After', 'FontSize', 18, 'FontWeight', 'bold');
% xlabel('ΔE_{00}', 'FontSize', labelFontSize, 'FontWeight', 'bold');
% ylabel('Pixel Count', 'FontSize', labelFontSize, 'FontWeight', 'bold');
% xlim([0 40]); grid on;
% set(gca, 'FontSize', tickFontSize, 'FontWeight', 'bold');
% 
% % Save histogram figure
% hist_out_path = fullfile(outputDir, [film_base_name '_deltaE_histograms.png']);
% % saveas(fig_hist, hist_out_path);
% exportgraphics(fig_hist, hist_out_path, 'Resolution', 300, 'BackgroundColor', 'none', 'ContentType', 'image');
% 
% fprintf('Saved histogram figure to: %s\n', hist_out_path);




%%
% Extract base name from film image path
[~, film_base_name, ~] = fileparts(film_path);

% Save corrected RGB image
rgb_out_path = fullfile(outputDir, [film_base_name '_corrected_rgb_srgb.tif']);
imwrite(corrected_rgb_srgb_img, rgb_out_path);
fprintf('Saved corrected RGB (sRGB) image to: %s\n', rgb_out_path);



%% Delta E stats

% Crop ΔE map interactively
disp('Select a rectangular region to crop in ΔE map (corrected vs. before ageing)...');
figure, imagesc(deltaE_before_img), axis image off;
title('Select crop region for ΔE_{00} map');
crop_rect = round(getrect);  % [x, y, width, height]
close;

% Apply cropping
x = crop_rect(1); y = crop_rect(2);
w_crop = crop_rect(3); h_crop = crop_rect(4);

% Ensure within image bounds
x_end = min(x + w_crop - 1, size(deltaE_before_img, 2));
y_end = min(y + h_crop - 1, size(deltaE_before_img, 1));

cropped_de = deltaE_before_img(y:y_end, x:x_end);

% Compute stats on cropped area
fprintf('CROPPED ΔE2000 stats (corrected vs. BEFORE ageing):\n');
fprintf('Mean ΔE2000:   %.3f\n', mean(cropped_de(:)));
fprintf('Median ΔE2000: %.3f\n', median(cropped_de(:)));
fprintf('Max ΔE2000:    %.3f\n', max(cropped_de(:)));




%%


function X_rootpoly = root_poly_features(input_data)

    a = input_data(:,1);  % R or channel 1
    b = input_data(:,2);  % G or channel 2
    c = input_data(:,3);  % B or channel 3

    % Avoid negative roots and division by zero
    a = max(a, 0); b = max(b, 0); c = max(c, 0);

    X_rootpoly = [ ...
        ones(size(a)), ...
        a, b, c, ...
        sqrt(a), sqrt(b), sqrt(c), ...
        a.*b, a.*c, b.*c, ...
        sqrt(a.*b), sqrt(a.*c), sqrt(b.*c), ...
        a.^2, b.^2, c.^2, ...
        sqrt(abs(a.^2 + b.^2 + c.^2))];
  
end


