clear; close all;

% --- INPUT FILE ---
inputFile = '/home/oem/eliza/data/reflectance/registered/yoda_reflectance_before.mat';      
outFolder = '/home/oem/eliza/data/xyz_lab_rgb/hyspex';      

% --- Create output folder if needed ---
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% --- Load input .mat file ---
S = load(inputFile);
data_cube = flip(S.data_cube, 2);      % Horizontal flip
spec_mask1 = flip(S.spec_mask, 2);      % Flip mask the same way
wl        = S.wl;
%%
% --- Binning parameters ---
binSize = 2;
[H, W, B] = size(data_cube);
h = floor(H / binSize);
w = floor(W / binSize);
binned_cube = zeros(h, w, B, 'like', data_cube);
spec_mask = false(h, w);

for b = 1:B
    band = data_cube(:,:,b);
    band = band(1:h*binSize, 1:w*binSize);  % Crop
    band = reshape(band, binSize, h, binSize, w);
    band = permute(band, [2 4 1 3]);
    binned_cube(:,:,b) = mean(mean(band, 3), 4);
end

% --- Process mask with same binning ---
mask = spec_mask1(1:h*binSize, 1:w*binSize);
mask = reshape(mask, binSize, h, binSize, w);
mask = permute(mask, [2 4 1 3]);
spec_mask = mean(mean(mask, 3), 4) > 0.5;  % Thresholded average
%%
% --- Load spectral data ---
CMF_path = '../../../data/CIE2degCMFs_full.csv';
D50_path = '../../../data/CIE_D50.txt';
fullCMFs = importdata(CMF_path);
D50_SPD  = importdata(D50_path);

cmf_interp = interp1(fullCMFs(:,1), fullCMFs(:,2:4), wl, 'linear', 'extrap');
ill_interp = interp1(D50_SPD(:,1), D50_SPD(:,2), wl, 'linear', 'extrap');

% --- Normalize so perfect white gives Y = 100 ---
S = cmf_interp .* ill_interp;     % [B x 3]
k = 100 / sum(S(:,2));            % normalization constant

% --- Compute XYZ, Lab, RGB ---
refl = reshape(binned_cube, [], B);       % [N × B]
XYZ  = k * (refl * S);                    % [N × 3]
Lab  = xyz2lab(XYZ);                      % [N × 3]
RGB  = xyz2prophoto(XYZ ./ 100, true);    % gamma-corrected ProPhoto RGB
RGB_lin = xyz2prophoto(XYZ ./ 100, false);% linear RGB

% --- Reshape to image format ---
XYZ_img     = reshape(XYZ, h, w, 3);
Lab_img     = reshape(Lab, h, w, 3);
RGB_img     = reshape(RGB, h, w, 3);
RGB_lin_img = reshape(RGB_lin, h, w, 3);

% --- Show RGB preview ---
figure; imshow(mat2gray(RGB_img)); title('ProPhoto RGB Preview');

% --- Save output (images + binned mask) ---
[~, baseName, ~] = fileparts(inputFile);
outputFile = fullfile(outFolder, [baseName '_xyz.mat']);
save(outputFile, 'XYZ_img', 'Lab_img', 'RGB_img', 'RGB_lin_img', 'spec_mask', 'wl');
fprintf('Saved converted data to:\n%s\n', outputFile);

%%
% --- Load saved result and check RGB image ---
loaded = load(outputFile);

figure;
imshow(mat2gray(loaded.RGB_img));
title('Loaded RGB Image (Verification)');
