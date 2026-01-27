%% ================= USER CONFIG =================

% Set comparison mode
% true  = compare AFTER vs FILM
% false = compare AFTER vs BEFORE
% true  = AFTER vs FILM
% false = AFTER vs BEFORE
compareAfterVsFilm = true;

% -------- Paths --------
afterMat  = "/Volumes/Study/Thesis/data/captures/reflectance_reg/mat/yoda_reflectance_after_reg.mat";
beforeMat = "/Volumes/Study/Thesis/data/captures/reflectance_reg/mat/cactus_reflectance_before.mat";

filmMat   = "/Volumes/Study/Thesis/data/captures/film_registered/yoda_halogen_fuji_exp0_balanced_registered.mat";

binSize = 2;


% Downsampling / binning parameters (as you already had)
spatial_bin = 2;     % e.g. 2x2 spatial binning
spectral_bin = 1;    % spectral binning if used

%% ================= LOAD DATA =================

% Always load AFTER
[after_cube, wl_after] = load_hsi_cube(path_after);

if compareAfterVsFilm
    [ref_cube, wl_ref] = load_hsi_cube(path_film);
    disp('Mode: AFTER vs FILM');
else
    [ref_cube, wl_ref] = load_hsi_cube(path_before);
    disp('Mode: AFTER vs BEFORE');
end

%% ================= ALIGN WAVELENGTHS =================

% Interpolate reference cube to AFTER wavelengths
% (avoids huge interp arrays if done band-by-band)
ref_cube_interp = zeros(size(after_cube), 'like', after_cube);

for i = 1:size(after_cube,1)
    for j = 1:size(after_cube,2)
        ref_spec = squeeze(ref_cube(i,j,:));
        ref_cube_interp(i,j,:) = interp1( ...
            wl_ref, ref_spec, wl_after, 'linear', 'extrap');
    end
end

%% ================= SPATIAL BINNING =================

after_cube_bin = spatial_bin_cube(after_cube, spatial_bin);
ref_cube_bin   = spatial_bin_cube(ref_cube_interp, spatial_bin);

%% ================= SPECTRAL BINNING (OPTIONAL) =================

if spectral_bin > 1
    after_cube_bin = spectral_bin_cube(after_cube_bin, spectral_bin);
    ref_cube_bin   = spectral_bin_cube(ref_cube_bin, spectral_bin);
end

%% ================= FINAL VARIABLES USED BY METRICS =================

% These two cubes are what ALL metrics should use
cube_A = after_cube_bin;
cube_B = ref_cube_bin;


epsVal = 1e-12;

[h,w,B] = size(cube_A);

X = reshape(cube_A, [], B);
Y = reshape(cube_B, [], B);

valid = all(isfinite(X),2) & all(isfinite(Y),2);
N = size(X,1);

%% ---------- RMSE ----------
rmse_pix = nan(N,1);
rmse_pix(valid) = sqrt(mean((X(valid,:) - Y(valid,:)).^2, 2));
RMSE_map = reshape(rmse_pix, h, w);

%% ---------- MAE ----------
mae_pix = nan(N,1);
mae_pix(valid) = mean(abs(X(valid,:) - Y(valid,:)), 2);
MAE_map = reshape(mae_pix, h, w);

%% ---------- GFC (cosine similarity) ----------
gfc_pix = nan(N,1);
num = sum(X(valid,:) .* Y(valid,:), 2);
den = sqrt(sum(X(valid,:).^2,2)) .* sqrt(sum(Y(valid,:).^2,2)) + epsVal;
gfc_pix(valid) = num ./ den;
GFC_map = reshape(gfc_pix, h, w);

%% ---------- CGFC ----------
CGFC_map = 1 - GFC_map;

%% ---------- SAM ----------
sam_pix = nan(N,1);
cosang = max(min(gfc_pix(valid),1), -1);
sam_pix(valid) = acosd(cosang);
SAM_map = reshape(sam_pix, h, w);

%% ---------- SCM (mean-centered correlation) ----------
scm_pix = nan(N,1);
for k = find(valid)'
    x = X(k,:);  y = Y(k,:);
    x0 = x - mean(x);
    y0 = y - mean(y);

    scm_pix(k) = (x0 * y0') / ...
        (sqrt(sum(x0.^2)) * sqrt(sum(y0.^2)) + epsVal);
end
SCM_map = reshape(scm_pix, h, w);

%% ---------- SID ----------
sid_pix = nan(N,1);
for k = find(valid)'
    p = X(k,:);  q = Y(k,:);
    p = p / (sum(p) + epsVal);
    q = q / (sum(q) + epsVal);

    sid_pix(k) = ...
        sum(p .* log((p + epsVal) ./ (q + epsVal))) + ...
        sum(q .* log((q + epsVal) ./ (p + epsVal)));
end
SID_map = reshape(sid_pix, h, w);

%% ---------- EDCS ----------
edcs_pix = nan(N,1);
for k = find(valid)'
    x = X(k,:);  y = Y(k,:);
    x = x / (sum(x) + epsVal);
    y = y / (sum(y) + epsVal);

    edcs_pix(k) = norm(cumsum(x) - cumsum(y));
end
EDCS_map = reshape(edcs_pix, h, w);
EDCS_norm = EDCS_map / max(EDCS_map(:));


figure('Name','2x4 Spectral Metrics','Color','w');
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');

nexttile;
imagesc(RMSE_map); axis image off; colorbar;
title('RMSE');

nexttile;
imagesc(MAE_map); axis image off; colorbar;
title('MAE');

nexttile;
imagesc(SAM_map,[0 30]); axis image off; colorbar;
title('SAM (deg)');

nexttile;
imagesc(CGFC_map,[0 0.3]); axis image off; colorbar;
title('CGFC');

nexttile;
imagesc(GFC_map,[0 1]); axis image off; colorbar;
title('GFC');

nexttile;
imagesc(SID_map,[0 prctile(SID_map(valid),99)]); axis image off; colorbar;
title('SID');

nexttile;
imagesc(EDCS_norm,[0 1]); axis image off; colorbar;
title('EDCS (norm)');

nexttile;
imagesc(SCM_map,[-1 1]); axis image off; colorbar;
title('SCM');




function cube_bin = spatial_bin_cube(cube, bin)
    if bin == 1
        cube_bin = cube;
        return;
    end

    [H,W,B] = size(cube);
    H2 = floor(H/bin);
    W2 = floor(W/bin);

    cube_bin = zeros(H2, W2, B, 'like', cube);

    for b = 1:B
        tmp = reshape(cube(1:H2*bin, 1:W2*bin, b), ...
                      bin, H2, bin, W2);
        tmp = squeeze(mean(mean(tmp,1),3));
        cube_bin(:,:,b) = tmp;
    end
end


function cube_bin = spectral_bin_cube(cube, bin)
    if bin == 1
        cube_bin = cube;
        return;
    end

    [H,W,B] = size(cube);
    B2 = floor(B/bin);
    cube_bin = zeros(H,W,B2,'like',cube);

    for k = 1:B2
        idx = (k-1)*bin + (1:bin);
        cube_bin(:,:,k) = mean(cube(:,:,idx), 3);
    end
end
