

clc; clear; close all;

%% ---------------- Load Hyperspectral Data ----------------
cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_120f-velvia-f8.hdr";
refFile  = "../../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";
rng(10);

hcube = hypercube(cubeFile);
inCUBE = hcube.DataCube;
bands = hcube.Wavelength;
[m, n, bd] = size(inCUBE);

hcube_ref = hypercube(refFile);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;
[mr, nr, bd_ref] = size(refCUBE);

% Flatten
lincube = reshape(inCUBE, [], bd);
lincube_ref = reshape(refCUBE, [], bd_ref);

%% ---------------- Compute XYZ from Spectral Data ----------------
% Load illuminant and CMFs
ill = importdata("../../../data/CIE_D65.txt");
CMFs_1931 = importdata("../../../data/CIE2degCMFs_1931.txt");
CMFs = CMFs_1931;

% Interpolate illuminant & CMFs to match captured wavelengths
illIP = interp1(ill(:,1), ill(:,2), bands, 'spline');
CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,3), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,4), bands, 'spline')];
sp_tristREF = CMFsIP .* illIP;

% Compute XYZ for input
xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2),1);

% Do the same for reference
illIP_ref = interp1(ill(:,1), ill(:,2), bands_ref, 'spline');
CMFsIP_ref = [interp1(CMFs(:,1), CMFs(:,2), bands_ref, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,3), bands_ref, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,4), bands_ref, 'spline')];
sp_tristREF_ref = CMFsIP_ref .* illIP_ref;

xyz_ref = (lincube_ref * sp_tristREF_ref) ./ sum(sp_tristREF_ref(:,2),1);

%% ---------------- Reshape to 10x14 (plus border) for patch indexing ----------------
% The chart is said to be 10×14 with a 1 patch grayscale border => effectively 12×16 total
% We'll reshape back so we can index patches by row, col
% (Adapt to your actual arrangement if different.)
chart_rows = 10;  % 1 patch border top & bottom
chart_cols = 14;  % 1 patch border left & right
xyz_input_2d = reshape(xyz_input, m, n, 3);  
xyz_ref_2d   = reshape(xyz_ref,   mr, nr, 3);

%% ---------------- Step 1: LUT/Tone Correction Using Grayscale Patches ----------------
% 1) Identify grayscale patches = border + center region (rows 5..6, cols 5..10)
% We'll gather their average XYZ from input & reference
grayRows = [1, chart_rows];       % top & bottom border
grayCols = [1, chart_cols];       % left & right border

patchXYZ_in = [];
patchXYZ_ref = [];

% Collect border patches
for row = 1:chart_rows
    for col = 1:chart_cols
        if row==1 || row==chart_rows || col==1 || col==chart_cols
            patchVal_in = xyz_input_2d(row, col, :);
            patchVal_ref = xyz_ref_2d(row, col, :);
            patchXYZ_in  = [patchXYZ_in; reshape(patchVal_in,1,3)];
            patchXYZ_ref = [patchXYZ_ref; reshape(patchVal_ref,1,3)];
        end
    end
end

% Collect center region (rows 5..6, cols 5..10)
for row = 5:6
    for col = 5:10
        patchVal_in = xyz_input_2d(row, col, :);
        patchVal_ref = xyz_ref_2d(row, col, :);
        patchXYZ_in  = [patchXYZ_in; reshape(patchVal_in,1,3)];
        patchXYZ_ref = [patchXYZ_ref; reshape(patchVal_ref,1,3)];
    end
end

% Create a simple 1D LUT or tone curve in each channel based on X, Y, Z
% For example, we can do an interpolation in each channel
correctedXYZ_2d_step1 = xyz_input_2d;
for ch = 1:3
    inVals  = patchXYZ_in(:,ch);
    refVals = patchXYZ_ref(:,ch);
    % Sort for monotonic interpolation
    [inVals_s, idx] = sort(inVals);
    refVals_s = refVals(idx);

    % Build interpolation function
    lutFun = @(x) interp1(inVals_s, refVals_s, x, 'linear', 'extrap');

    % Apply to that channel in the entire image
    channelData = xyz_input_2d(:,:,ch);
    channelCorrected = lutFun(channelData);
    correctedXYZ_2d_step1(:,:,ch) = channelCorrected;
end

 correctedXYZ_step1 = reshape(correctedXYZ_2d_step1, [], 3);

%% ---------------- Step 2: Another Correction (3rd-degree polynomial) ----------------

% Flatten the reference as well
xyz_ref_flat = reshape(xyz_ref_2d, [], 3);

% (Optional) train-test split
num_samples = size(correctedXYZ_step1, 1);
perm = randperm(num_samples);
trainCount = round(0.8*num_samples);
train_idx = perm(1:trainCount);
test_idx  = perm(trainCount+1:end);

rgb_corr1 = xyz2rgb(correctedXYZ_step1, 'ColorSpace', 'prophoto-rgb');
rgb_ref   = xyz2rgb(xyz_ref,   'ColorSpace', 'prophoto-rgb');

rgb_input_train = rgb_corr1(train_idx, :);
rgb_ref_train   = rgb_ref(train_idx, :);


% Build 3rd-degree polynomial
X_poly_train = poly3_features(rgb_input_train);
coeffs_poly  = pinv(X_poly_train) * rgb_ref_train;

% Apply to full data
X_poly_full = poly3_features(rgb_corr1);
rgb_corr2 = X_poly_full * coeffs_poly;

%% ---------------- Evaluate & Display ----------------
% Convert to Lab for error evaluation
lab_ref = xyz2lab(xyz_ref_flat);
lab_step2 = rgb2lab(rgb_corr2, "ColorSpace", "prophoto-rgb");

% Evaluate
% Suppose you have an evaluate_error function
mRows = m;  % or chart_rows
nCols = n;  % or chart_cols
output_folder = '../../../results/error_maps';
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

evaluate_error(lab_ref, lab_step2, test_idx, mRows, nCols, '2step_', output_folder, "chart_2step_correction.png");


% Flatten back to apply next correction
rgb_corr1 = reshape(rgb_corr1, [], 3);
rgb_corr2 = reshape(rgb_corr1, [], 3);


% Plot comparison images side by side with correct aspect ratio (10x14 chart)
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Ensure the patches don't look squished
tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');

% Original XYZ
nexttile;
imagesc(reshape(xyz_input, 10, 14, 3));
axis equal off;
title('Original XYZ');

% After Step 1 (LUT correction)
nexttile;
imagesc(reshape(rgb_corr1, 10, 14, 3));
axis equal off;
title('After Step1 LUT');

% After Step 2 (Polynomial correction)
nexttile;
imagesc(reshape(rgb_corr2, 10, 14, 3));
axis equal off;
title('After Step2 Poly');

% Add row and column labels
for i = 1:3
    nexttile(i);
    set(gca, 'XTick', 1:14, 'XTickLabel', 1:14, ...
             'YTick', 1:10, 'YTickLabel', 1:10, ...
             'TickLength', [0 0], 'FontSize', 10);
end

% Adjust the layout for better visualization
sgtitle('Two-Step Color Correction: XYZ -> LUT -> Polynomial');


%% Polynomial Feature Function
function X_poly = poly3_features(input_data)
    % For a 3-column input [a, b, c], create 3rd-degree expansions
    a = input_data(:,1);
    b = input_data(:,2);
    c = input_data(:,3);

    feat1 = a; feat2 = b; feat3 = c;
    feat4 = a.^2; feat5 = b.^2; feat6 = c.^2;
    feat7 = a.*b; feat8 = a.*c; feat9 = b.*c;
    feat10 = a.^3; feat11 = b.^3; feat12 = c.^3;
    feat13 = a.^2 .* b; feat14 = a.^2 .* c; feat15 = b.^2 .* a;
    feat16 = b.^2 .* c; feat17 = c.^2 .* a; feat18 = c.^2 .* b;
    feat19 = a.*b.*c;

    X_poly = [feat1, feat2, feat3, feat4, feat5, feat6, feat7, feat8, feat9, ...
              feat10, feat11, feat12, feat13, feat14, feat15, feat16, feat17, feat18, feat19];
end
