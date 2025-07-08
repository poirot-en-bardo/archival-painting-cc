clear; close all;

% --- INPUT FILE ---
inputFile = '/home/oem/eliza/data/film_scans/registered/yoda_led_kodak_exp0_balanced_registered.mat';      
outFolder = '/home/oem/eliza/data/xyz_lab_rgb/film';      

% --- Create output folder if needed ---
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% --- Load registered cube ---
S = load(inputFile);
cube = flip(S.wb_cube_registered, 2);  % Horizontal flip
wl   = S.wl;
%%
% --- Binning ---
binSize = 2;
[H, W, B] = size(cube);
h = floor(H / binSize);
w = floor(W / binSize);
binned_cube = zeros(h, w, B, 'like', cube);

for b = 1:B
    band = cube(:,:,b);
    band = band(1:h*binSize, 1:w*binSize);
    reshaped = reshape(band, binSize, h, binSize, w);
    reshaped = permute(reshaped, [2 4 1 3]);
    binned_cube(:,:,b) = mean(mean(reshaped, 3), 4);
end

% --- Load spectral data ---
CMF_path = '../../../data/CIE2degCMFs_full.csv';
D50_path = '../../../data/CIE_D50.txt';
LED_path = '../../../data/film/CVCL10bands.txt';

fullCMFs = importdata(CMF_path);       % [λ, x y z]
D50_SPD  = importdata(D50_path);       % [λ, SPD]
LEDset   = readmatrix(LED_path, 'Delimiter','\t');
wl_led   = LEDset(:,1);                % λ axis for LED SPDs
led_spds = LEDset(:,2:end);            % each column is LED SPD
bandN    = size(led_spds,2);
%
% --- Integrate CMFs and D50 over LED bands ---
d50_spd  = interp1(D50_SPD(:,1), D50_SPD(:,2), wl_led, 'linear', 'extrap');
cmf_xyz  = interp1(fullCMFs(:,1), fullCMFs(:,2:4), wl_led, 'linear', 'extrap');

d50_band = zeros(bandN, 1);
cmf_band = zeros(bandN, 3);
for k = 1:bandN
    led_spd = led_spds(:,k);
    norm_led = led_spd / sum(led_spd);
    d50_band(k)   = sum(d50_spd .* norm_led);
    cmf_band(k,:) = sum(cmf_xyz .* norm_led, 1);
end

% --- Normalize so white gives Y = 100 ---
k = 100 / sum(d50_band .* cmf_band(:,2));

% --- Compute XYZ, Lab, RGB ---
radCube = double(binned_cube) .* reshape(d50_band, 1, 1, []);
linCube = reshape(radCube, [], bandN);     % [N × B]
XYZ     = k * (linCube * cmf_band);        % [N × 3]
Lab  = xyz2lab(XYZ);
RGB  = xyz2prophoto(XYZ./100, true);
RGB_lin = xyz2prophoto(XYZ ./ 100, false);

% --- Reshape to image format ---
XYZ_img     = reshape(XYZ, h, w, 3);
Lab_img     = reshape(Lab, h, w, 3);
RGB_img     = reshape(RGB, h, w, 3);
RGB_lin_img = reshape(RGB_lin, h, w, 3);

% --- Show RGB preview ---
figure; imshow(mat2gray(RGB_img)); title('ProPhoto RGB Preview');
%%
% --- Save output ---
[~, baseName, ~] = fileparts(inputFile);
outputFile = fullfile(outFolder, [baseName '_xyz.mat']);
save(outputFile, 'XYZ_img', 'Lab_img', 'RGB_img', 'RGB_lin_img', 'wl');
fprintf('Saved converted data to:\n%s\n', outputFile);
%%
% --- Optional: Verification plot ---
loaded = load(outputFile);

% Show ProPhoto RGB image
figure;
imshow(mat2gray(loaded.RGB_img));
title('Loaded ProPhoto RGB Image (Verification)');

% Compute sRGB from XYZ (normalize to [0,1] first)
XYZ_norm = loaded.XYZ_img ./ 100;
sRGB_img = xyz2rgb(XYZ_norm, 'ColorSpace', 'srgb', 'WhitePoint', 'd50');
sRGB_img = min(max(sRGB_img, 0), 1);  % Clamp

% Show sRGB image
figure;
imshow(sRGB_img);
title('Computed sRGB from XYZ (Verification)');
