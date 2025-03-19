
clc; clear; close all;
roof = double(intmax('uint16'));
histoFACT = 200; % Histogram normalization factor

%% Load spectral data
ill = importdata('../../../data/CIE_D65.txt'); % Illuminant
CMFs_1931 = importdata('../../../data/CIE2degCMFs_1931.txt');
CMFs_2006 = importdata('../../../data/CIE2degCMFs_2006.txt');
CMFs = CMFs_1931;
rng(10);

%% Select the reference and the cube to be corrected
cube_folder = '../../../data/colorChecker_SG/colorChecker_SG_Elias';
cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_fuji-frame4.hdr";
refFile  = "../../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

%% Load and process the cube to be corrected (input cube)
hcube = hypercube(cubeFile);
inCUBE = hcube.DataCube;
bands = hcube.Wavelength;
[m, n, bd] = size(inCUBE);
lincube = reshape(inCUBE, [], bd);

% Interpolate illuminant and CMFs to captured wavelengths
illIP = interp1(ill(:,1), ill(:,2), bands, 'spline');
CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,3), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,4), bands, 'spline')];
sp_tristREF = CMFsIP .* illIP;
% Compute XYZ values for the input cube (normalized by Y)
xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);

%%% Load and process the reference cube
hcube_ref = hypercube(refFile);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;
[m, n, bd] = size(refCUBE);
lincube_ref = reshape(refCUBE, [], bd);

illIP = interp1(ill(:,1), ill(:,2), bands_ref, 'spline');
CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands_ref, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,3), bands_ref, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,4), bands_ref, 'spline')];
sp_tristREF = CMFsIP .* illIP;
% Compute XYZ values for the reference cube
xyz_ref = (lincube_ref * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);

num_samples = size(xyz_input, 1);
perm = randperm(num_samples);
train_idx = perm(1:round(0.8*num_samples));
test_idx  = perm(round(0.8*num_samples)+1:end);

%% ----------------- CAM16-Based Regression -----------------
% Define viewing conditions for CAM16 conversion:
whiteXYZ = [95.047, 100.000, 108.883]; % D65 white point
LA = 20;   % Adapting luminance (cd/m^2)
F  = 1.0;  % Surround factor (average surround)

% Convert the input and reference XYZ values to CAM16 JMh values
cam16_input = xyz2cam16(xyz_input, whiteXYZ, LA, F);
cam16_ref   = xyz2cam16(xyz_ref, whiteXYZ, LA, F);

% For regression, split into training and testing data as before
cam16_input_train = cam16_input(train_idx, :);
cam16_ref_train   = cam16_ref(train_idx, :);

% Compute 3rd-degree polynomial features for CAM16 training data
X_poly_cam16_train = poly3_features(cam16_input_train);

% Compute regression coefficients using pseudoinverse
coeffs_cam16 = pinv(X_poly_cam16_train) * cam16_ref_train;

% Apply the regression to the full dataset
X_poly_cam16 = poly3_features(cam16_input);
corrected_cam16 = X_poly_cam16 * coeffs_cam16;

%% ----------------- Evaluate & Display Results -----------------
% Convert the corrected CAM16 values back to XYZ
corrected_xyz_from_cam16 = cam16_to_xyz(corrected_cam16, whiteXYZ, LA, F);

% Convert both the reference and corrected XYZ values to Lab for error evaluation
lab_from_xyz_ref  = xyz2lab(xyz_ref);
lab_from_xyz_cor  = xyz2lab(corrected_xyz_from_cam16);


output_folder = '../../../results/error_maps';
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
[~, img_name, ~] = fileparts(cubeFile);


% Evaluate and display the error maps (using your existing evaluate_error function)
evaluate_error(lab_from_xyz_ref, lab_from_xyz_cor, test_idx, m, n, 'CAM16_', output_folder, img_name + "_error_cam16.png");


%% Conversion functions
function cam16 = xyz2cam16(xyz, whiteXYZ, LA, F)
% xyz2cam16 - Convert from XYZ to CAM16 JMh values
%
% Syntax:
%   cam16 = xyz2cam16(xyz, whiteXYZ, LA, F)
%
% Inputs:
%   xyz      - n-by-3 matrix of XYZ values
%   whiteXYZ - 1-by-3 vector for the white point (e.g., D65)
%   LA       - Adapting field luminance in cd/m^2 (e.g., 20)
%   F        - Surround factor (typical values: 1.0 for average surround)
%
% Output:
%   cam16    - n-by-3 matrix containing [J, M, h] (Lightness, 
%              colorfulness, hue angle in degrees)
%
% Reference: Li, et al. “CAM16: A Color Appearance Model for the New ICC Profile”
%

% Constants for CAT16 conversion matrix
M_CAT16 = [ 0.401288,  0.650173, -0.051461;
           -0.250268,  1.204414,  0.045854;
           -0.002079,  0.048952,  0.953127];

% Convert the input white point to LMS
LMS_white = (M_CAT16 * whiteXYZ')';

% Convert XYZ to LMS for the input
LMS = (M_CAT16 * xyz')';

% Calculate the degree of adaptation D (using the formula from CAM16)
% Here we use a common approximation. Adjust as needed.
D = F * (1 - (1/3.6)*exp(-(LA + 42)/92));
D = min(max(D, 0), 1);  % constrain D to [0, 1]

% Compute the adapted LMS values
LMS_c = D .* (whiteXYZ(1) ./ LMS_white(1)) + (1-D);
M_c   = D .* (whiteXYZ(2) ./ LMS_white(2)) + (1-D);
S_c   = D .* (whiteXYZ(3) ./ LMS_white(3)) + (1-D);
LMS_a = LMS .* repmat([LMS_c, M_c, S_c], size(LMS,1), 1);

% Nonlinear response compression parameters
FL = (0.2 * (LA^0.4)) + 0.1; % Simplified; actual FL computation is more involved

% Apply non-linear response compression (simplified version)
F_L = FL;
LMS_prime = (400 * (LMS_a./100).^(0.42)) ./ (27.13 + (LMS_a./100).^(0.42)) + 0.1;

% Compute the achromatic response A for each channel (simplified)
A = sum(LMS_prime, 2);

% Compute J (lightness correlate)
% Here, we use a simplified formula. The full CAM16 model uses a more complex equation.
J = 100 * (A / sum(LMS_white)).^0.69;

% Compute chroma correlate (M) and hue angle (h)
% For demonstration, we compute M as a simple function of the differences between channels.
% A full implementation would use opponent color dimensions.
a = LMS_prime(:,1) - LMS_prime(:,2);
b = LMS_prime(:,2) - LMS_prime(:,3);
M = sqrt(a.^2 + b.^2);
h = atan2d(b, a);
h(h < 0) = h(h < 0) + 360;  % Ensure hue is in [0, 360]

% Return CAM16 correlates: [J, M, h]
cam16 = [J, M, h];
end


function xyz = cam16_to_xyz(cam16, whiteXYZ, LA, F)
% cam16_to_xyz - Approximate inverse conversion from CAM16 JMh values to XYZ
%
% Syntax:
%   xyz = cam16_to_xyz(cam16, whiteXYZ, LA, F)
%
% Inputs:
%   cam16    - n-by-3 matrix containing [J, M, h] (Lightness, colorfulness, hue angle in degrees)
%   whiteXYZ - 1-by-3 vector for the white point (e.g., D65: [95.047, 100.000, 108.883])
%   LA       - Adapting field luminance in cd/m^2 (e.g., 20)
%   F        - Surround factor (typical value: 1.0 for average surround)
%
% Output:
%   xyz      - n-by-3 matrix of approximated XYZ values
%
% Note: This implementation is a rough inversion of the simplified xyz2cam16 function
% provided earlier and is intended for demonstration purposes only.

% Unpack CAM16 correlates
J = cam16(:,1); % Lightness
M = cam16(:,2); % Colorfulness
h = cam16(:,3); % Hue angle (degrees)

% Constants and CAT16 matrix (as used in the forward conversion)
M_CAT16 = [ 0.401288,  0.650173, -0.051461;
           -0.250268,  1.204414,  0.045854;
           -0.002079,  0.048952,  0.953127];
       
% Convert the white point to LMS
LMS_white = (M_CAT16 * whiteXYZ')';

% --- Step 1: Recover the achromatic response A ---
% From the forward function, J was computed as:
%   J = 100 * (A / sum(LMS_white)).^0.69
% Invert to recover A:
A = (J/100).^(1/0.69) * sum(LMS_white);

% --- Step 2: Estimate opponent dimensions (a and b) ---
% In the forward function, a and b were approximated as differences between
% the non-linear responses of the channels. Here we assume that M (the colorfulness)
% distributes into opponent dimensions based on the hue angle:
a = M .* cosd(h);
b = M .* sind(h);

% --- Step 3: Approximate the non-linear compressed responses ---
% We need to approximate LMS_prime such that their sum approximates A.
% Here we assume a simple split:
%   LMS_prime(:,1) = A/3 + a/2
%   LMS_prime(:,2) = A/3
%   LMS_prime(:,3) = A/3 - b/2
LMS_prime = zeros(size(cam16,1), 3);
LMS_prime(:,1) = A/3 + a/2;
LMS_prime(:,2) = A/3;
LMS_prime(:,3) = A/3 - b/2;

% --- Step 4: Invert the non-linear response compression ---
% The forward compression was:
%   LMS_prime = 400 * (LMS_a/100).^0.42 ./ (27.13 + (LMS_a/100).^0.42) + 0.1;
% We solve for LMS_a (approximate inversion):
temp = (LMS_prime - 0.1) ./ (400 - (LMS_prime - 0.1)); 
LMS_a = 100 * (27.13 * temp).^(1/0.42);

% --- Step 5: Undo the adaptation scaling ---
% In the forward function, adaptation was applied as:
%   LMS_a = LMS .* repmat([LMS_c, M_c, S_c], size(LMS,1), 1)
% where for each channel:
%   channel_factor = D * (whiteXYZ(channel)/LMS_white(channel)) + (1-D)
% First, compute D as in the forward function:
D = F * (1 - (1/3.6)*exp(-(LA+42)/92));
D = min(max(D, 0), 1);
LMS_c = D * (whiteXYZ(1) / LMS_white(1)) + (1-D);
M_c   = D * (whiteXYZ(2) / LMS_white(2)) + (1-D);
S_c   = D * (whiteXYZ(3) / LMS_white(3)) + (1-D);

% Undo the scaling:
LMS = zeros(size(LMS_a));
LMS(:,1) = LMS_a(:,1) / LMS_c;
LMS(:,2) = LMS_a(:,2) / M_c;
LMS(:,3) = LMS_a(:,3) / S_c;

% --- Step 6: Invert the CAT16 transformation to get back to XYZ ---
% Solve for XYZ: xyz = inv(M_CAT16) * LMS'
xyz = (M_CAT16 \ LMS')';
end

function X_poly = poly3_features(input_data)
    % Constructs a 3rd-degree polynomial expansion (without constant term)
    % for a 3-column input. For input columns a, b, c, this returns:
    % [ a, b, c, a.^2, b.^2, c.^2, a.*b, a.*c, b.*c, ...
    %   a.^3, b.^3, c.^3, a.^2.*b, a.^2.*c, b.^2.*a, b.^2.*c, c.^2.*a, c.^2.*b, a.*b.*c ]
    a = input_data(:,1);
    b = input_data(:,2);
    c = input_data(:,3);
    
    % Degree 1
    feat1 = a;
    feat2 = b;
    feat3 = c;
    
    % Degree 2
    feat4 = a.^2;
    feat5 = b.^2;
    feat6 = c.^2;
    feat7 = a.*b;
    feat8 = a.*c;
    feat9 = b.*c;
    
    % Degree 3
    feat10 = a.^3;
    feat11 = b.^3;
    feat12 = c.^3;
    feat13 = a.^2 .* b;
    feat14 = a.^2 .* c;
    feat15 = b.^2 .* a;
    feat16 = b.^2 .* c;
    feat17 = c.^2 .* a;
    feat18 = c.^2 .* b;
    feat19 = a .* b .* c;
    
    X_poly = [feat1, feat2, feat3, feat4, feat5, feat6, feat7, feat8, feat9, ...
              feat10, feat11, feat12, feat13, feat14, feat15, feat16, feat17, feat18, feat19];
end