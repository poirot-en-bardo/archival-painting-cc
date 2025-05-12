hdr_path = "/home/oem/eliza/data/processed/yoda_halogen_reflectance_full.hdr";

cube_data = hypercube(hdr_path);
wl  = cube_data.Wavelength;               % assume same for both
cube = cube_data.DataCube;

[rows, cols, B] = size(cube);
linCube = double( reshape(cube, rows*cols, B) );  

%% Load and interpolate illuminant + CMFs ———
ill = importdata('../../../data/CIE_D65.txt');          
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

%% Convert to RGB —————


%sRGB
sRGB = xyz2rgb(XYZimg, 'ColorSpace','srgb');


%% ProPhoto-RGB 
pPhoto = xyz2rgb(XYZimg, 'ColorSpace','prophoto-rgb');


%% 3) Scale to uint16 and save
% sRGB → uint16 (0→0, 1→65535)
sRGB16 = uint16(sRGB * 65535);

% ProPhoto → uint16
pPhoto = uint16(pPhoto * 65535);

%%
imwrite(sRGB16, 'plots/yoda_halogen_sRGB.png');
imwrite(pPhoto, 'plots/yoda_halogen_pPhoto.png');

%%
figure; imshow(sRGB16);
title('sRGB (16-bit)');
%%
figure; imshow(pPhoto);
title('ProPhoto (16-bit)');