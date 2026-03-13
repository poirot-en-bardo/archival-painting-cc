%% ============================================================
%  Paths
%% ============================================================
close; clear all;
% %%
% afterMat = "/Volumes/Study/Thesis/data/captures/reflectance_reg/mat/cactus_reflectance_after_reg.mat";
% % afterMat = "/Volumes/Study/Thesis/data/captures/reflectance_reg/mat/yoda_reflectance_after_reg.mat";
% 
% beforeMat = "/Volumes/Study/Thesis/data/captures/reflectance_reg/mat/cactus_reflectance_before.mat";
% % beforeMat = "/Volumes/Study/Thesis/data/captures/reflectance_reg/mat/yoda_reflectance_before.mat";
% % filmMat  = "/Volumes/Study/Thesis/data/captures/film_registered/cactus_halogen_kodak_exp0_balanced_registered.mat";
% % filmMat  = "/Volumes/Study/Thesis/data/captures/film_registered/cactus_led_fuji_exp0_balanced_registered.mat";
% % filmMat  = "/Volumes/Study/Thesis/data/captures/film_registered/cactus_led_fuji_underexp_balanced_registered.mat";
% % filmMat  = "/Volumes/Study/Thesis/data/captures/film_registered/yoda_halogen_fuji_exp0_balanced_registered.mat";
% % filmMat  = "/Volumes/Study/Thesis/data/captures/film_registered/yoda_led_kodak_exp0_balanced_registered.mat";
% filmMat = "/Volumes/Study/Thesis/data/captures/film_registered/yoda_halogen_fuji_underexp_balanced_registered.mat";
% 




%% ============================================================
%  Load after_ageing hypercube
%% ============================================================
% Ca = load(afterMat);
% Cb = load(beforeMat);
% spec_mask_after = Ca.spec_mask;
% spec_mask_before = Cb.spec_mask;
% 
% binSize = 2;




%%

% aa_cube = Ca.data_cube;
% aa_wl = Ca.wl;

% bb_cube = Cb.data_cube;
% bb_wl = Cb.wl;



%% ============================================================
%  Limit spectral range (380–780 nm)
%% ============================================================
% wlMask = aa_wl >= 380 & aa_wl <= 780;
% 
% aa_cube = aa_cube(:,:,wlMask);
% aa_wl   = aa_wl(wlMask);




% wlMask = bb_wl >= 380 & bb_wl <= 780; 
% bb_cube = bb_cube(:,:,wlMask);
% bb_wl   = bb_wl(wlMask);

% %% ============================================================
% %  Load film (10-band reference)
% %% ============================================================
% S = load(filmMat);
% 
% film = S.wb_cube_registered;   % H × W × 10
% wl   = S.wl(:);                % 10 wavelengths
% 
% %%
% % % Indices to keep
% % keepIdx = [3 6 9];
% % 
% % % Subselect film cube
% % film = film(:,:,keepIdx);   % now H × W × 3
% % 
% % % Subselect wavelengths
% % wl   = wl(keepIdx);         % now 3 × 1
% 
% 
% %% ============================================================
% %  Spectral downsampling: after_ageing → film wavelengths
% %% ============================================================
% [H, W, ~] = size(aa_cube);
% 
% aa_rs = reshape(aa_cube, [], numel(aa_wl));   % [N × B]
% 
% aa_down = interp1( ...
%     aa_wl, ...
%     aa_rs.', ...
%     wl, ...
%     'linear', 'extrap' ...
% ).';
% 
% aa_down = reshape(aa_down, H, W, numel(wl));  % H × W × 10
% 
% %%

% 
% [H, W, ~] = size(bb_cube);
% bb_rs = reshape(bb_cube, [], numel(bb_wl));   % [N × B]
% 
% bb_down = interp1( ...
%     bb_wl, ...
%     bb_rs.', ...
%     wl, ...
%     'linear', 'extrap' ...
% ).';
% 
% bb_down = reshape(bb_down, H, W, numel(wl));  % H × W × 10
% 
% %% ============================================================
% %  Spatial binning (same logic as your code)
% %% ============================================================
% after_bin   = bin_hypercube(aa_down, binSize);
% film_bin = bin_hypercube(film,    binSize);
% before_bin= bin_hypercube(bb_down,    binSize);
%%


afterMat = "/Volumes/Study/Thesis/data/captures/binned/cactus_after_bin.mat";
% afterMat = "/Volumes/Study/Thesis/data/captures/binned/yoda_after_bin.mat";
beforeMat = "/Volumes/Study/Thesis/data/captures/binned/cactus_before_bin.mat";
% beforeMat = "/Volumes/Study/Thesis/data/captures/binned/yoda_before_bin.mat";
% filmMat = "/Volumes/Study/Thesis/data/captures/binned/cactus_halogen_kodak_exp0_bin.mat";
% filmMat = "/Volumes/Study/Thesis/data/captures/binned/cactus_led_fuji_underexp_bin.mat";
filmMat = "/Volumes/Study/Thesis/data/captures/binned/cactus_led_fuji_exp0_bin.mat";
% filmMat = "/Volumes/Study/Thesis/data/captures/binned/yoda_halogen_fuji_overexp_bin.mat";
% filmMat = "/Volumes/Study/Thesis/data/captures/binned/yoda_led_kodak_exp0_bin.mat";
% filmMat = "/Volumes/Study/Thesis/data/captures/binned/yoda_halogen_fuji_exp0_bin.mat";




%%
struct_after = load(afterMat);
struct_before = load(beforeMat);
after_bin = struct_after.after_cube_bin;
before_bin = struct_before.before_cube_bin;
struct_film = load(filmMat);
film_bin = struct_film.film_cube_bin;
wl = struct_film.wl;
spec_mask_after = struct_after.spec_mask;
spec_mask_before = struct_before.spec_mask;


binSize = 2;
compareAfterVsFilm = false;
epsVal = 1e-12;


%%



if compareAfterVsFilm == true
    compare_bin = film_bin;
else
    compare_bin = before_bin;
end



%% ============================================================
[h, w, B] = size(after_bin);   % B = 10

X = reshape(after_bin,   [], B);   % [N × 10]
Y = reshape(compare_bin, [], B);   % [N × 10]

valid = all(isfinite(X),2) & all(isfinite(Y),2);
% 
rmse_pix = nan(size(X,1),1);
rmse_pix(valid) = sqrt(mean((X(valid,:) - Y(valid,:)).^2, 2));

RMSE_map = reshape(rmse_pix, h, w);

mae_pix = nan(size(X,1),1);
mae_pix(valid) = mean(abs(X(valid,:) - Y(valid,:)), 2);

MAE_map = reshape(mae_pix, h, w);



% 
% Standard GFC (cosine similarity)
gfc_pix = nan(size(X,1),1);
num = sum(X(valid,:) .* Y(valid,:), 2);
den = sqrt(sum(X(valid,:).^2,2)) .* sqrt(sum(Y(valid,:).^2,2));
gfc_pix(valid) = num ./ den;
GFC_map = reshape(gfc_pix, h, w);

% Complemented GFC (CGFC) = 1 - GFC
CGFC_map = 1 - GFC_map;


% SAM
sam_pix = nan(size(X,1),1);
sam_pix(valid) = acosd( max(min(gfc_pix(valid),1), -1) );  % <--- use gfc_pix, not cgfc_pix
SAM_map = reshape(sam_pix, h, w);
%%

%% ============================================================
%  Hybrid RMSE + SAM (50/50)
%% ============================================================
if compareAfterVsFilm == true
    mask = spec_mask_after;
else
    mask = spec_mask_before | spec_mask_after;
end
mask = logical(mask);
RMSE_map(mask) = 0;
SAM_map(mask) = 0;

% Robust normalization (avoid outliers dominating)
percent = 95;

rmse_scale = prctile(RMSE_map(:), percent);
sam_scale  = prctile(SAM_map(:),  percent);
cgfc_scale = prctile(CGFC_map(:),  percent);

RMSE_norm = RMSE_map / rmse_scale;
SAM_norm  = SAM_map  / sam_scale;
CGFC_norm = CGFC_map / cgfc_scale;

% RMSE_norm = RMSE_map / max(RMSE_map(:));
% SAM_norm  = SAM_map  / max(SAM_map(:));
% CGFC_norm = CGFC_map / cgfc_scale;
% 
% SAM_unit = SAM_map / 90;
% SAM_norm = min( SAM_unit / prctile(SAM_unit(:),95), 1 );



% Clip to [0,1] for stability
RMSE_norm = min(max(RMSE_norm, 0), 1);
SAM_norm  = min(max(SAM_norm,  0), 1);
CGFC_norm  = min(max(CGFC_norm,  0), 1);



%  Multiplicative Hybrid (AND-like)
%
HYBRID_MULT = RMSE_norm .* SAM_norm;

HYBRID_MULT(mask) = 0;

% figure('Name','Hybrid Multiplicative','Color','w');
% imagesc(HYBRID_MULT);
% axis image off; colorbar; colormap("parula");
% title('Hybrid Multiplicative (RMSE × SAM)');


%% ===== Hybrid Multiplicative (publication style) =====

% -------- Styling parameters --------
figSize   = [0.1, 0.2, 1, 0.8];   % [left, bottom, width, height]
fontSize  = 28;
labelSize = 28;

% Optional: set color limits for consistency across figures
% clim_range = [0 1];   % uncomment if you want fixed scaling

% -------- Create figure --------
fig = figure('Name','Hybrid Multiplicative', ...
             'Units','normalized', ...
             'Position', figSize, ...
             'Color','w');

imagesc(HYBRID_MULT);
axis image off;
colormap("parula");

% If you want fixed limits:
% clim(clim_range);

% -------- Colorbar (big & bold) --------
cb = colorbar;
cb.FontSize   = labelSize;
cb.FontWeight = 'bold';
% cb.LineWidth  = 2;

% cb.Label.String     = 'Hybrid (RMSE × SAM)';
cb.Label.FontSize   = labelSize;
cb.Label.FontWeight = 'bold';
clim([0 1]);                 % make sure data range is 0..1
cb.Ticks = 0:0.2:1;
cb.TickLabels = string(0:0.2:1);   % optional, for clean formatting

% -------- Axes font (for consistency) --------
set(gca, 'FontSize', fontSize, 'FontWeight', 'bold');

% title('Hybrid Multiplicative', 'FontSize', fontSize, 'FontWeight', 'bold');

%% ===== Save figure (tight, with colorbar, no padding) =====

outFolder = "/Volumes/Study/Thesis/thesis-repo/results/plots/change_detection/spectral_hybrid/new/multiplicative";

if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

[~, filmName, ~] = fileparts(filmMat);
filmName = string(filmName);

if contains(filmName, "_bin")
    baseName = extractBefore(filmName, "_bin");
else
    baseName = filmName;
end

if compareAfterVsFilm == false
    baseName = baseName + "_GT";
else
    baseName = baseName + "_hybrid";
end


outFile = fullfile(outFolder, baseName + ".jpg");   % PNG better for papers

% exportgraphics(fig, outFile, ...
%     'Resolution', 300, ...
%     'BackgroundColor', 'none', ...
%     'ContentType', 'image');

fprintf('Saved Hybrid figure to:\n%s\n', outFile);




%%
% 
% figure;
% imagesc(SAM_map, [0 30]);
% axis image off;
% colorbar; colormap("parula");
% % title('SAM');
% cb = colorbar;
% 
% % ---- Make colorbar text big & bold ----
% cb.FontSize = 30;        % super big
% cb.FontWeight = 'bold'; % bold numbers



% 
% figure;
% imagesc(RMSE_map, [0 0.2]);
% axis image off;
% colorbar; colormap("parula");
% % title('RMSE');
% cb = colorbar;
% 
% % ---- Make colorbar text big & bold ----
% cb.FontSize = 30;        % super big
% cb.FontWeight = 'bold'; % bold numbers



%%

% -------- User-defined output folder --------
outFolder = "/Volumes/Study/Thesis/thesis-repo/results/plots/change_detection/spectral_hybrid/new/metrics";   % <-- change if needed

% Make sure folder exists
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% -------- Build base filename from filmMat --------
[~, filmName, ~] = fileparts(filmMat);
filmName = string(filmName);

% Remove '_bin' and everything after
if contains(filmName, "_bin")
    baseName = extractBefore(filmName, "_bin");
else
    baseName = filmName;
end

% Append _GT if comparing against GT (before)
baseName = baseName + "_SAM";


% Full output file
outFile = fullfile(outFolder, baseName + ".jpg");

% -------- Tight export (image + colorbar, no white padding) --------
% exportgraphics(gcf, outFile, ...
%     'Resolution', 300, ...
%     'BackgroundColor', 'none');

fprintf('Saved figure to:\n%s\n', outFile);


%%

% figure;
% imagesc(RMSE_norm, [0 1]);
% axis image off;
% colorbar; colormap("parula");
% title('RMSE per pixel');
% 
% figure;
% imagesc(SAM_norm, [0 1]);
% axis image off;
% colorbar; colormap("parula");
% title('SAM per pixel');
% 
% % 
% % 
% % 
% % 

% % 
% figure;
% imagesc(CGFC_map, [0 1]);
% axis image off;
% colorbar;
% title('CGFC');



% 
% % 50/50 weighted hybrid
% w_rmse = 0.5;
% w_sam  = 0.5;
% 
% HYBRID_RMSE_SAM = ...
%     w_rmse * RMSE_norm + ...
%     w_sam  * SAM_norm;
% 

% Visualization
% figure('Name','Hybrid RMSE + SAM','Color','w');
% imagesc(HYBRID_RMSE_SAM, [0 1]);
% axis image off; colorbar; colormap("parula");
% title('Hybrid Metric (0.5 RMSE + 0.5 SAM)');
% 
% w_rmse = 0.5;
% w_cgfc  = 0.5;
% 
% HYBRID_RMSE_CGFC = ...
%     w_rmse * RMSE_norm + ...
%     w_sam  * CGFC_norm;
% 
% % Visualization
% figure('Name','Hybrid RMSE + SAM','Color','w');
% imagesc(HYBRID_RMSE_CGFC, [0 1]);
% axis image off; colorbar; colormap("parula");
% title('Hybrid Metric (0.5 RMSE + 0.5 CGFC)');





%% ============================================================
%  Helper: spatial binning (identical to your implementation)
%% ============================================================
function binned_cube = bin_hypercube(data_cube, binSize)

    [H, W, B] = size(data_cube);
    h = floor(H / binSize);
    w = floor(W / binSize);

    binned_cube = zeros(h, w, B, 'like', data_cube);

    for b = 1:B
        band = data_cube(:,:,b);
        band = band(1:h*binSize, 1:w*binSize);

        band = reshape(band, binSize, h, binSize, w);
        band = permute(band, [2 4 1 3]);
        binned_cube(:,:,b) = mean(mean(band, 3), 4);
    end
end

function binned_img = bin_image_2d(img, binSize)

    [H, W] = size(img);
    h = floor(H / binSize);
    w = floor(W / binSize);

    img = img(1:h*binSize, 1:w*binSize);

    img = reshape(img, binSize, h, binSize, w);
    img = permute(img, [2 4 1 3]);

    binned_img = mean(mean(img, 3), 4);
end


%%






%%
%% ============================================================
%  Interactive ROI spectral comparison (HSI vs Film)
%% ============================================================





% 
% bandToShow = 5;   % <-- choose band for visualization (1–10)
% 
% figure('Name','Select ROI (double-click to confirm)','Color','w');
% imagesc(after_bin(:,:,bandToShow));
% axis image off;
% colorbar;
% title(sprintf('HSI band %d – zoom/pan, then draw ROI (double-click)', bandToShow));
% 
% % Enable zoom + pan before drawing
% zoom on;
% pan on;
% uiwait(msgbox('Zoom/pan as needed, then click OK to draw ROI','ROI','modal'));
% zoom off;
% pan off;
% 
% % Draw ROI (double-click inside rectangle to confirm)
% roi = drawrectangle('Color','r');
% wait(roi);
% 
% mask = createMask(roi);
% 
% %% ============================================================
% %  Extract and average spectra inside ROI
% %% ============================================================
% 
% % HSI spectra
% roi_hsi = reshape(after_bin, [], B);
% roi_film = reshape(compare_bin, [], B);
% 
% mask_vec = mask(:);
% 
% mean_hsi  = mean(roi_hsi(mask_vec,:), 1, 'omitnan');
% mean_film = mean(roi_film(mask_vec,:), 1, 'omitnan');
% 
% %% ============================================================
% %  Plot spectra
% %% ============================================================
% % 
% figure('Name','ROI Mean Spectra','Color','w');
% plot(wl, mean_hsi, '-o','LineWidth',2); hold on;
% plot(wl, mean_film,'-s','LineWidth',2);
% grid on;
% 
% xlabel('Wavelength (nm)');
% ylabel('Mean reflectance');
% legend('HSI (after ageing)','HSI (before ageing)','Location','best');
% title(sprintf('Mean spectrum)', bandToShow));


%%



% Flatten SAM map to a vector
sam_vals = SAM_map(:);

% Remove NaNs
sam_vals = sam_vals(isfinite(sam_vals));

% Plot histogram
% figure('Color','w');
% histogram(sam_vals, 50, 'Normalization', 'probability'); % 50 bins, normalized to probability
% xlabel('SAM angle (degrees)');
% ylabel('Probability');
% title('Histogram of SAM values');
% grid on;


% Flatten CGFC map to a vector
cgfc_vals = CGFC_map(:);

% Remove NaNs
cgfc_vals = cgfc_vals(isfinite(cgfc_vals));

% % Plot histogram
% figure('Color','w');
% histogram(cgfc_vals, 50, 'Normalization', 'probability');
% xlabel('CGFC (1 - cosine similarity)');
% ylabel('Probability');
% title('Histogram of CGFC values');
% grid on;
% 
% %%
% cgfc_vals = SAM_map(:);
% 
% % Histogram
% figure;
% histogram(cgfc_vals, 100);
% xlabel('SAM');
% ylabel('Number of pixels');
% title('SAM Distribution');
% 
% cgfc_vals = RMSE_map(:);
% 
% % Histogram
% figure;
% histogram(cgfc_vals, 100);
% xlabel('RMSE');
% ylabel('Number of pixels');
% title('RMSE Distribution');

% % CDF
% figure;
% cdfplot(cgfc_vals);
% xlabel('SAM');
% ylabel('Cumulative probability');
% title('CDF of SAM values');

%%
% %% ============================================================
%visulaise spectra

bandToShow = 5;   % <-- choose band for visualization (1–10)

% Display the selected band from AFTER cube for ROI selection
figure('Name','Select ROI (double-click to confirm)','Color','w');
imagesc(after_bin(:,:,bandToShow));
axis image off;
colorbar;
title(sprintf('Select ROI on AFTER cube, Band %d (double-click)', bandToShow));

% Enable zoom + pan before drawing ROI
zoom on; pan on;
uiwait(msgbox('Zoom/pan as needed, then click OK to draw ROI','ROI','modal'));
zoom off; pan off;

% Draw ROI
roi = drawrectangle('Color','r');
wait(roi);

mask = createMask(roi);
mask_vec = mask(:);

%% ============================================================
% Extract and average spectra inside ROI
%% ============================================================

% Reshape cubes
[H,W,B] = size(after_bin);
roi_after  = reshape(after_bin, [], B);
roi_before = reshape(before_bin, [], B);
roi_film   = reshape(film_bin, [], B);

% Compute mean spectra inside ROI
mean_after  = mean(roi_after(mask_vec,:), 1, 'omitnan');
mean_before = mean(roi_before(mask_vec,:), 1, 'omitnan');
mean_film   = mean(roi_film(mask_vec,:), 1, 'omitnan');

%% ============================================================
% Plot all three spectra together
%% ============================================================
% figure('Name','ROI Mean Spectra','Color','w'); hold on;
% 
% plot(wl, mean_after,  '-o', 'LineWidth',2, 'DisplayName','After ageing');
% plot(wl, mean_before, '-s', 'LineWidth',2, 'DisplayName','Before ageing');
% % plot(wl, mean_film,   '-d', 'LineWidth',2, 'DisplayName','Film reference');
% 
% grid on;
% xlabel('Wavelength (nm)');
% ylabel('Mean reflectance');
% legend('Location','best');
% % title(sprintf('Mean spectra', bandToShow));


%% ============================================================
% Plot all three spectra together (THESIS STYLE)
%% ============================================================

fig = figure('Name','ROI Mean Spectra', ...
             'Color','w', ...
             'Units','normalized', ...
             'Position',[0.2 0.2 0.6 0.6]);
hold on;

% ---- Super thick, bold lines & big markers ----
p1 = plot(wl, mean_after,  '-s', ...
    'Color', [0.8500 0.3250 0.0980], ...
    'LineWidth', 6, ...
    'MarkerSize', 10, ...
    'DisplayName','After ageing');

p2 = plot(wl, mean_before, '-s', ...
    'Color', [0 0.4470 0.7410], ...
    'LineWidth', 6, ...
    'MarkerSize', 10, ...
    'DisplayName','Before ageing');

% Optional film reference
% p3 = plot(wl, mean_film,   '-d', ...
%     'LineWidth', 4, ...
%     'MarkerSize', 10, ...
%     'DisplayName','Film reference');

% ---- Axes & grid styling ----
ylim([0.02 0.07]);

ax = gca;
ax.FontSize   = 33;
ax.FontWeight = 'bold';
ax.LineWidth  = 3;
ax.TickLength = [0.015 0.015];

grid on;
ax.GridLineWidth = 2;
ax.MinorGridLineWidth = 1.5;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';

% ---- Labels & legend ----
xlabel('Wavelength (nm)', 'FontSize', 33, 'FontWeight', 'bold');
ylabel('Mean reflectance', 'FontSize', 33, 'FontWeight', 'bold');

leg = legend('Location','best');
leg.FontSize   = 33;
leg.FontWeight = 'bold';
leg.Box = 'off';

% title('Mean spectra (ROI)', 'FontSize', 30, 'FontWeight', 'bold');

%% ===== Save plot (based on afterMat, ultramarine) =====

outFolder = "/Volumes/Study/Thesis/thesis-repo/results/plots/roi_spectra";

if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% -------- Build base filename from afterMat --------
[~, afterName, ~] = fileparts(afterMat);
afterName = string(afterName);

% Extract before '_after'
if contains(afterName, "_after")
    baseName = extractBefore(afterName, "_after");
else
    baseName = afterName;  % fallback
end

% Append required suffix
baseName = baseName + "_red";

outFile = fullfile(outFolder, baseName + ".png");

exportgraphics(fig, outFile, ...
    'Resolution', 300, ...
    'BackgroundColor', 'none', ...
    'ContentType', 'image');

fprintf('Saved ROI spectra plot to:\n%s\n', outFile);
