clc; clear; close all;
rng(10);

%% ---------------- Load Hyperspectral Data ----------------
cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_120f-velvia-f8.hdr";
refFile  = "../../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

hcube = hypercube(cubeFile);
inCUBE = hcube.DataCube;
bands = hcube.Wavelength;
[m, n, bd] = size(inCUBE);

hcube_ref = hypercube(refFile);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;
[mr, nr, bd_ref] = size(refCUBE);

% Flatten spectral data
lincube = reshape(inCUBE, [], bd);
lincube_ref = reshape(refCUBE, [], bd_ref);

%% ---------------- Compute XYZ from Spectral Data ----------------
% Load illuminant and CMFs
ill = importdata("../../../data/CIE_D65.txt");
CMFs_1931 = importdata("../../../data/CIE2degCMFs_1931.txt");
CMFs = CMFs_1931;

% For input cube:
illIP = interp1(ill(:,1), ill(:,2), bands, 'spline');
CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,3), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,4), bands, 'spline')];
sp_tristREF = CMFsIP .* illIP;
xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2),1);

% For reference cube:
illIP_ref = interp1(ill(:,1), ill(:,2), bands_ref, 'spline');
CMFsIP_ref = [interp1(CMFs(:,1), CMFs(:,2), bands_ref, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,3), bands_ref, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,4), bands_ref, 'spline')];
sp_tristREF_ref = CMFsIP_ref .* illIP_ref;
xyz_ref = (lincube_ref * sp_tristREF_ref) ./ sum(sp_tristREF_ref(:,2),1);

%% ---------------- Reshape to 10x14 ----------------
% The chart is arranged as 10 rows x 14 columns.
xyz_input_2d = reshape(xyz_input, 10, 14, 3);
xyz_ref_2d   = reshape(xyz_ref,   10, 14, 3);

%% ---------------- Convert XYZ to RGB ----------------
% Convert using the 'prophoto-rgb' color space
rgb_input = xyz2rgb(xyz_input, 'ColorSpace', 'prophoto-rgb');
rgb_ref   = xyz2rgb(xyz_ref,   'ColorSpace', 'prophoto-rgb');

% Reshape to 10x14 for display and patch extraction
rgb_input_2d = reshape(rgb_input, 10, 14, 3);
rgb_ref_2d   = reshape(rgb_ref, 10, 14, 3);

%% ---------------- Extract White and Black Patches ----------------
% White patch at (1,1) and Black patch at (2,1)
input_white = squeeze(rgb_input_2d(1,1,:))';    % 1x3 vector
input_black = squeeze(rgb_input_2d(2,1,:))';

ref_white = squeeze(rgb_ref_2d(1,1,:))';
ref_black = squeeze(rgb_ref_2d(2,1,:))';

%% ---------------- Build Tone Curve LUT ----------------
% For each channel, create a simple linear mapping from the input tone curve
% using two anchor points: black and white.
% The function maps the input value (in the range [input_black, input_white])
% to the reference value (in the range [ref_black, ref_white]).
% We'll do this per channel.
toneCurve = cell(1,3);  % cell array to store LUT functions for R, G, and B.
for ch = 1:3
    inVals = [input_black(ch), input_white(ch)];   % anchor points from input
    refVals = [ref_black(ch), ref_white(ch)];        % corresponding target values
    toneCurve{ch} = @(x) interp1(inVals, refVals, x, 'linear', 'extrap');
end

%% ---------------- Apply Tone Curve to Entire Image ----------------
% Flatten the RGB image for processing.
rgb_input_flat = reshape(rgb_input_2d, [], 3);
rgb_corrected_flat = zeros(size(rgb_input_flat));
for i = 1:size(rgb_input_flat,1)
    for ch = 1:3
        rgb_corrected_flat(i,ch) = toneCurve{ch}(rgb_input_flat(i,ch));
    end
end

% Reshape the corrected flat data back to 10x14 image
rgb_corrected_2d = reshape(rgb_corrected_flat, 10, 14, 3);

% Scale the values from [0, 1] to [0, 65535] and convert to uint16
maxVal = double(intmax('uint16'));  % This is 65535 for uint16
rgb_corrected_uint16 = uint16(rgb_corrected_2d * maxVal);

% Save the corrected image as a TIFF file (using no compression)
% imwrite(rgb_corrected_uint16, 'corrected_image_uint16_prophoto.tif', 'Compression', 'none');


%% Display Results
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
imagesc(rgb_input_2d);
axis equal off;
title('Original Input (RGB)');

nexttile;
imagesc(rgb_ref_2d);
axis equal off;
title('Reference (RGB)');

nexttile;
imagesc(reshape(rgb_corrected_flat, 10, 14, 3));
axis equal off;
title('Corrected with Tone Curve (RGB)');

sgtitle('Single-Step Tone Curve Correction Based on White & Black Anchors');

%% Visualize the Tone Curve in 3D (Optional)
% For each channel, we can visualize the 1D LUT.
% figure;
% for ch = 1:3
%     subplot(1,3,ch);
%     x = linspace(0, 1, 100);
%     y = toneCurve{ch}(x);
%     plot(x, y, 'LineWidth', 2);
%     xlabel('Input Value'); ylabel('Corrected Value');
%     title(sprintf('Channel %d', ch));
%     grid on;
% end
% sgtitle('1D Tone Curves for Each RGB Channel');

