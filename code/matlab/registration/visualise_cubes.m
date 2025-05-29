hdr_path = "/Volumes/School/Thesis/data/Baby Yoda/Yoda_reflectance.hdr";

cube_data = hypercube(hdr_path);
wl  = cube_data.Wavelength;               % assume same for both
cube = cube_data.DataCube;

% Find indices within 380–780 nm
valid_idx = find(wl >= 380 & wl <= 780);

% Subset cube and wavelengths
cube = cube(:, :, valid_idx);
wl = wl(valid_idx);

[rows, cols, B] = size(cube);
linCube = double( reshape(cube, rows*cols, B) );  

%% Load and interpolate illuminant + CMFs ———
ill = importdata('../../../data/CIE_D50.txt');          
CMFs = importdata('../../../data/CIE2degCMFs_1931.txt'); 

% Interpolate to your band centers
illIP = interp1(ill(:,1),     ill(:,2), wl, 'spline');          % [B×1]
CMFx  = interp1(CMFs(:,1), CMFs(:,2), wl, 'spline');           % [B×1]
CMFy  = interp1(CMFs(:,1), CMFs(:,3), wl, 'spline');
CMFz  = interp1(CMFs(:,1), CMFs(:,4), wl, 'spline');
CMFsIP = [CMFx, CMFy, CMFz];                                    % [B×3]

% Pre-multiply by illuminant to get spectral tristimulus functions
S = CMFsIP .* illIP;    % [B×3]

%% Compute XYZ
k = sum(S(:,2), 1);               % normalizing constant
XYZ = (linCube * S) ./ k;        % [N×3]

% Reshape back to image
XYZimg = reshape(XYZ, rows, cols, 3);


%% chromatic adaptation
% 
% % %% --- Chromatic-adapt XYZ (D65 ➜ D50) with the Bradford CAT -------------
% wpSrc = whitepoint('D65');         % [0.95047 1.00000 1.08883]
% wpDst = whitepoint('D50');         % [0.96422 1.00000 0.82521]
% 
% M_Bradford = [ 0.8951  0.2664 -0.1614 ;
%               -0.7502  1.7135  0.0367 ;
%                0.0389 -0.0685  1.0296 ];
% 
% % build 3×3 CAT
% LMS_src =  M_Bradford * wpSrc';
% LMS_dst =  M_Bradford * wpDst';
% D       =  diag(LMS_dst ./ LMS_src);
% M_CAT   =  inv(M_Bradford) * D * M_Bradford;     % D65 → D50
% 
% % --- reshape → multiply → reshape back ----------------------------------
% XYZ_flat      = reshape(XYZimg, [], 3);   % [N×3]   (N = rows*cols)
% XYZ_D50_flat  = XYZ_flat * M_CAT.';       % post-multiply each row
% XYZ_D50       = reshape(XYZ_D50_flat, size(XYZimg));
% pPhoto = xyz2rgb(XYZ_D50, 'ColorSpace','prophoto-rgb', 'OutputType','uint16');

%% Convert to RGB —————
%sRGB
sRGB = xyz2rgb(XYZimg, 'ColorSpace','srgb', 'OutputType', 'uint16', 'WhitePoint','d50');

% ProPhoto-RGB 
pPhoto = xyz2rgb(XYZimg, 'ColorSpace','prophoto-rgb', 'OutputType', 'uint16','WhitePoint','d50');

%linear RGB
linear = xyz2rgb(XYZimg, 'ColorSpace', 'linear-rgb', 'OutputType', 'uint16', 'WhitePoint','d50');


%% 3) Scale to uint16 and save
% sRGB → uint16 (0→0, 1→65535)
% sRGB16 = uint16(sRGB * 65535);
% 
% % ProPhoto → uint16
% pPhoto = uint16(pPhoto * 65535);

%%
outDir  = fullfile('/home/oem/eliza/mac-shared/HSI/before');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
outFile = fullfile(outDir, 'yoda_D50_sRGB.png');
% imwrite(sRGB, outFile);
%%
outFile = fullfile(outDir, 'yoda_D50_pPhoto.png');
imwrite(pPhoto, outFile);
%%
outFile = fullfile(outDir, 'yoda_D50_linear.png');
% imwrite(linear, outFile);

%%
figure; imshow(sRGB);
title('sRGB (16-bit)');
%%
figure; imshow(pPhoto);
title('ProPhoto (16-bit)');

%%
figure; imshow(linear);
title('Linear RGB (16-bit)');