clear; close all;

% --- Filepaths ---
matPath   = "/home/oem/eliza/data/film_scans/white_balanced/yoda_led_kodak_exp0_balanced.mat";
hdrFix    = "/home/oem/eliza/data/processed/reflectance/registered/yoda_reflectance_full_before.hdr";
refBandHS = 100;
refBandMS = 10;

% --- Load cubes ---
S       = load(matPath, 'wb_cube', 'wl');
movCube = flip(S.wb_cube, 1);             % [rows × cols × bands]
wlMS    = S.wl;

hcFix   = hypercube(hdrFix);
fixCube = gather(hcFix.DataCube);         % [rows × cols × bands]
wlFix   = hcFix.Wavelength;
%%
% --- Select ROI on MOVING cube (MSI) ---
threshold = 1;
bandImg = double(movCube(:,:,refBandMS));
bandImg(bandImg > threshold) = threshold;
bandImg = mat2gray(bandImg);
bandImgAdj = imadjust(bandImg);

figure; imshow(bandImgAdj); title('Select ROI on MOVING MSI band');
hRec = drawrectangle('StripeColor','r'); wait(hRec);
pos = round(hRec.Position);
close;

movCrop = movCube(pos(2):(pos(2)+pos(4)-1), pos(1):(pos(1)+pos(3)-1), :);
fixCrop = fixCube;
%%
% --- Prepare registration bands ---
movRef = double(movCrop(:,:,refBandMS));
fixRef = double(fixCrop(:,:,refBandHS));

movRef = mat2gray(min(movRef, threshold));
fixRef = mat2gray(min(fixRef, threshold));
hp = imsharpen(movRef, 'Radius', 1, 'Amount', 3);
movRefAdj = imadjust(hp, [], [], 0.5);
fixRefAdj = imadjust(fixRef);
%%
% --- Load or select control points
movingPointsFile = 'movingPoints.mat';
fixedPointsFile  = 'fixedPoints.mat';

if isfile(movingPointsFile) && isfile(fixedPointsFile)
    load(movingPointsFile, 'movingPoints');
    load(fixedPointsFile,  'fixedPoints');
    [movingPoints, fixedPoints] = cpselect(movRefAdj, fixRefAdj, movingPoints, fixedPoints, 'Wait', true);
else
    [movingPoints, fixedPoints] = cpselect(movRefAdj, fixRefAdj, 'Wait', true);
end

save(movingPointsFile, 'movingPoints');
save(fixedPointsFile,  'fixedPoints');
%%
% --- Fit geometric transform
transformType = 'polynomial';  
tform = fitgeotrans(movingPoints, fixedPoints, transformType,4);

% --- Apply to all bands
outputView = imref2d(size(fixRefAdj));
[~, ~, nBands] = size(movCrop);
registeredCube = zeros([size(fixRefAdj), nBands], 'like', movCrop);

for b = 1:nBands
    registeredCube(:,:,b) = imwarp(movCrop(:,:,b), tform, 'OutputView', outputView);
end
%
% --- Show overlay
regBand = registeredCube(:,:,refBandMS);
figure; imshowpair(fixRefAdj, regBand);
title('Overlay: Fixed vs Registered (Manual cpselect)');

%% save
% --- Generate output filename based on MSI matPath
[~, baseName, ~] = fileparts(matPath);
baseName = char(baseName);  % Ensure it's a character vector
outputFolder = '/home/oem/eliza/data/film_scans/registered';

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

outputFile = fullfile(outputFolder, [baseName '_registered.mat']);

% Save
wb_cube_registered = registeredCube;
wl = wlMS;
save(outputFile, 'wb_cube_registered', 'wl', '-v7.3');
fprintf('Registered cube saved to:\n%s\n', outputFile);

%% check saved file
% --- Reload the saved registered cube to confirm
loaded = load(outputFile);  % Loads wb_cube_registered and wl
regCube = loaded.wb_cube_registered;

% --- Extract the registered band and original fixed band
regBandFinal = mat2gray(min(regCube(:,:,refBandMS), 1));
fixBandFinal = mat2gray(min(fixCrop(:,:,refBandHS), 1));  % already cropped

% --- Contrast enhancement (optional)
regBandAdj = imadjust(regBandFinal, [], [], 0.6);
fixBandAdj = imadjust(fixBandFinal);
%
% --- Show blended overlay
figure;
imshowpair(fixBandAdj, regBandAdj);
title('Overlay from Reloaded Registered Cube vs. Fixed Band');


%%
%% --- Compute D50 sRGB renderings and display montage ---

% --- Load spectral data only once
D50_path      = '../../../data/CIE_D50.txt';
CMF_path      = '../../../data/CIE2degCMFs_1931.txt';
CMF_path_full = '../../../data/CIE2degCMFs_full.csv';
led_path      = '../../../data/film/CVCL10bands.txt';

D50      = importdata(D50_path);               % [λ, SPD]
CMFs     = importdata(CMF_path);               % For HSI
CMFs_f   = importdata(CMF_path_full);          % For MSI

LEDset   = readmatrix(led_path, 'Delimiter','\t');  
led_wl   = LEDset(:,1);             
led_spds = LEDset(:,2:end);         
bandN    = size(led_spds,2);

% --- HSI (fixed) cube → D50 RGB ---
illIP = interp1(D50(:,1), D50(:,2), wlFix, 'spline');
CMFx  = interp1(CMFs(:,1), CMFs(:,2), wlFix, 'spline');
CMFy  = interp1(CMFs(:,1), CMFs(:,3), wlFix, 'spline');
CMFz  = interp1(CMFs(:,1), CMFs(:,4), wlFix, 'spline');
S_hsi = [CMFx, CMFy, CMFz] .* illIP;
k_hsi = sum(S_hsi(:,2));
linFix = reshape(fixCrop, [], size(fixCrop,3));
XYZ_hsi = (linFix * S_hsi) / k_hsi;
XYZimg_hsi = reshape(XYZ_hsi, size(fixCrop,1), size(fixCrop,2), 3);
RGB_hsi = xyz2rgb(XYZimg_hsi, 'ColorSpace','srgb', 'WhitePoint','d50');
RGB_hsi = min(max(RGB_hsi, 0), 1);

% --- Registered MSI cube → D50 RGB ---
d50_spd = interp1(D50(:,1), D50(:,2), led_wl, 'linear', 'extrap');
cmf_xyz = interp1(CMFs_f(:,1), CMFs_f(:,2:4), led_wl, 'linear', 'extrap');
d50_band = zeros(bandN,1);
cmf_band = zeros(bandN,3);
for k = 1:bandN
    norm_spd = led_spds(:,k) / sum(led_spds(:,k));
    d50_band(k)   = sum(d50_spd .* norm_spd);
    cmf_band(k,:) = sum(cmf_xyz .* norm_spd,1);
end
data_k_norm = sum(d50_band .* cmf_band(:,2));
radMSI  = registeredCube .* reshape(d50_band, 1, 1, bandN);
linMSI  = reshape(radMSI, [], bandN) * cmf_band;
XYZ_msi = linMSI / data_k_norm;
XYZimg_msi = reshape(XYZ_msi, size(registeredCube,1), size(registeredCube,2), 3);
RGB_msi = xyz2rgb(XYZimg_msi, 'ColorSpace','srgb', 'WhitePoint','d50');
RGB_msi = min(max(RGB_msi, 0), 1);

% --- Display RGBs in montage
figure;
imshowpair(RGB_hsi, RGB_msi, 'montage');
title('D50 sRGB Rendering — Left: Fixed HSI, Right: Registered MSI', 'FontSize', 14);
