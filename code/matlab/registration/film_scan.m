ms_img_folder = "/home/oem/eliza/data/film_scans/before_ageing/normal/yoda_led_kodak_exp0/FlatFielded";
led_path = '../../../data/film/CVCL10bands.txt';
cmf_path = '../../../data/CIE2degCMFs_full.csv';

roof = double(intmax('uint16'));
bandN = 10;


% Get all TIFF filenames in the folder:
files = dir(fullfile(ms_img_folder, '*.tif'));

% Read the first image to grab size + class:
firstImg = imread(fullfile(ms_img_folder, files(1).name));
[H, W, C] = size(firstImg);


% Preallocate cube (H×W×7) of the same class:
cube = zeros(H, W, numel(files), class(firstImg));

% Loop through and fill:
for k = 1:numel(files)
    img = imread(fullfile(ms_img_folder, files(k).name));
    cube(:,:,k) = img;
end

%% rendering the rgb image
% imshow(cube(:,:,3), []);

D50_SPD = importdata('../../../data/CIE_D50.txt');  

LEDset = readmatrix(led_path, 'Delimiter','\t');
wl_led = LEDset(:,1);       

fullCMFs = importdata(cmf_path);
%% downsampling CMFs to LED set
wl_cmf = fullCMFs(:, 1);
idx = wl_cmf >= 380 & wl_cmf <= 780;
wl_cmf = wl_cmf(idx);      
cmfs_vis     = fullCMFs(idx, :);  
CMFs = (cmfs_vis(:,2:end)'*LEDset(:,2:end)./sum(LEDset(:,2:end)))';


[r, c, ~] = size(cube);
lincube_capture = reshape(double(cube)/roof,[],bandN);
xyz_MS_full= (lincube_capture * CMFs)./sum(CMFs(:,2),1);

rgbImg_MS_full = xyz2rgb(xyz_MS_full,'ColorSpace','prophoto-rgb');
rgbImg_MS_full  = reshape(rgbImg_MS_full, r, c, 3);  

%% downsampling CMFs and D50 to LED set

bandN = size(LEDset, 2) - 1;    
wl_led = LEDset(:,1);           
d50_spd = interp1(D50_SPD(:,1), D50_SPD(:,2), wl_led, 'linear', 'extrap');
cmf_xyz = interp1(fullCMFs(:,1), fullCMFs(:,2:4), wl_led, 'linear', 'extrap');

d50_band = zeros(bandN, 1);    
cmf_band = zeros(bandN, 3);     

for k = 1:bandN
    led_spd = LEDset(:, k+1);    
    norm_led = led_spd / sum(led_spd);  
    
    d50_band(k)      = sum(d50_spd   .* norm_led);
    cmf_band(k, 1)   = sum(cmf_xyz(:,1) .* norm_led); 
    cmf_band(k, 2)   = sum(cmf_xyz(:,2) .* norm_led); 
    cmf_band(k, 3)   = sum(cmf_xyz(:,3) .* norm_led); 
end

cube_dbl = double(cube) / double(intmax('uint16'));  
cube_illum = cube_dbl .* reshape(d50_band, 1, 1, bandN);

[nRows, nCols, ~] = size(cube_illum);
cube_flat = reshape(cube_illum, [], bandN);

XYZ = cube_flat * cmf_band;         

k = sum(d50_band .* cmf_band(:,2));
XYZ = XYZ / k;

RGB = xyz2rgb(XYZ, 'ColorSpace', 'prophoto-rgb', 'WhitePoint','d50');
rgbImg_MS_full = reshape(RGB, nRows, nCols, 3);

imshow(rgbImg_MS_full);
title('ProPhoto RGB - Downsampled via LED SPDs');


%% D50 & CMFs interpolated to LED peaks
wl_led = LEDset(:,1);               
[~,~,bandN] = size(cube);    

wl_peaks = zeros(bandN,1);
for k = 1:bandN
    [~, idx] = max(LEDset(:,k+1));  % first col = wavelength
    wl_peaks(k) = wl_led(idx);      
end

% interpolate D50 and CMFs to the peak wavelengths
D50      = importdata('../../../data/CIE_D50.txt');
cmf      = importdata('../../../data/CIE2degCMFs_full.csv');
d50_IP   = interp1(D50(:,1), D50(:,2), wl_peaks, 'pchip', 'extrap');    
cmf_x    = interp1(cmf(:,1), cmf(:,2), wl_peaks, 'pchip', 'extrap');
cmf_y    = interp1(cmf(:,1), cmf(:,3), wl_peaks, 'pchip', 'extrap');
cmf_z    = interp1(cmf(:,1), cmf(:,4), wl_peaks, 'pchip', 'extrap');

CMFs_interp = [cmf_x, cmf_y, cmf_z];

S = CMFs_interp .* d50_IP;

% Normalising constant 
k = sum(d50_IP .* cmf_y);

% computing XYZ
cube_dbl = double(cube) / roof;  % [0–1]
[H, W, B] = size(cube_dbl);         
linCube  = reshape(cube_dbl, [], B);     
XYZ      = (linCube * S) / k;             
XYZimg   = reshape(XYZ, H, W, 3);
rgbImg_MS_full = xyz2rgb(XYZimg, 'ColorSpace', 'prophoto-rgb', 'WhitePoint','d50');
sRGB = xyz2rgb(XYZimg, 'ColorSpace', 'srgb', 'WhitePoint','d50');
sRGB = max(min(sRGB,1),0);
figure; imshow(sRGB); title('with D50');

%% checking illuminant effect --> CMFs + coloured illuminant
fake_illum = ones(size(wl_peaks));
fake_illum(wl_peaks < 550) = 0;  % reddish illuminant

S_fake = CMFs_interp .* fake_illum(:);
k = sum(fake_illum .* cmf_y);
XYZ_fake = (linCube * S_fake) / k;
XYZimg   = reshape(XYZ_fake, H, W, 3);

% illuminant white point
white_spectrum = ones(size(fake_illum));  
XYZ_white = (white_spectrum(:)' * S_fake) / k;  


% rgbImg_MS_full = xyz2rgb(XYZimg, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', XYZ_white);
% sRGB = xyz2rgb(XYZimg, 'ColorSpace', 'srgb', 'WhitePoint', XYZ_white); 

% rgbImg_MS_full = xyz2rgb(XYZimg, 'ColorSpace', 'prophoto-rgb', 'WhitePoint','d65');
sRGB = xyz2rgb(reshape(XYZ_fake, H, W, 3), 'ColorSpace', 'srgb'); % should give the white point of the illuminant actually
sRGB = max(min(sRGB,1),0);
figure; imshow(sRGB);
title('RGB with Blue Bands Zeroed');

%% only CMFs interpolated to LED peaks
cmf_x    = interp1(cmf(:,1), cmf(:,2), wl_peaks, 'pchip', 'extrap');
cmf_y    = interp1(cmf(:,1), cmf(:,3), wl_peaks, 'pchip', 'extrap');
cmf_z    = interp1(cmf(:,1), cmf(:,4), wl_peaks, 'pchip', 'extrap');
CMFs_interp = [cmf_x, cmf_y, cmf_z];

cube_dbl = double(cube) / roof;
[H, W, B] = size(cube_dbl);
linCube = reshape(cube_dbl, [], B);
S = CMFs_interp;
k = sum(cmf_y);

XYZ = (linCube * S) / k;
XYZimg = reshape(XYZ, H, W, 3);
sRGB = xyz2rgb(XYZimg, 'ColorSpace', 'srgb', 'WhitePoint','d50');
sRGB = max(min(sRGB,1),0);
figure;imshow(sRGB); title('just CMFs');



%% White balancing based on CC white patch
figure; imshow(rgbImg_MS_full);
title('Draw rectangle around white patch and double-click');
hW = drawrectangle();            
wait(hW);                 % block until double-click
whitePos = hW.Position;

xw    = max(1, floor(whitePos(1)));
yw    = max(1, floor(whitePos(2)));
ww    = floor(whitePos(3));
hw    = floor(whitePos(4));
xw_end = min(size(rgbImg_MS_full,2), xw + ww  - 1);
yw_end = min(size(rgbImg_MS_full,1), yw + hw  - 1);

patch_rgb   = rgbImg_MS_full(yw:yw_end, xw:xw_end, :); 
white_patch = squeeze( mean(mean(patch_rgb,1),2) );    % average

% white-balance the whole image
rgbWB = rgbImg_MS_full ./ reshape(white_patch,1,1,3);
rgbWB = max(min(rgbWB,1),0);   % [0,1]
rgb16 = im2uint16(rgbWB);



%% cropping the image
figure; imshow(rgb16);
title('Draw rectangle to crop and double-click');
hC = drawrectangle();            
wait(hC);                 
cropPos = hC.Position;    

% clamp & round
xc    = max(1, floor(cropPos(1)));
yc    = max(1, floor(cropPos(2)));
wc    = floor(cropPos(3));
hc2   = floor(cropPos(4));
xc_end = min(size(rgb16,2), xc + wc - 1);
yc_end = min(size(rgb16,1), yc + hc2 - 1);

rgbCropped = rgb16(yc:yc_end, xc:xc_end, :);


figure; imshow(rgbCropped);
title('Cropped RGB Image');
%% Save the cropped image
outDir  = fullfile('/home/oem/eliza/mac-shared/film_scans/before');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
outFile = fullfile(outDir, 'yoda_led_kodak_exp0_d65.png');

imwrite(rgbCropped, outFile);


%% Graphs for comparison

% CMFs interpolated to LED set peaks vs downsampled (scaled by the average)
figure;
subplot(1,3,1);
plot(wl_peaks, cmf_band(:,1), 'bo-', 'DisplayName','SPD-weighted X'); hold on;
plot(wl_peaks, cmf_x,        'rx--', 'DisplayName','Peak-interp X');
xlabel('Wavelength (nm)'); ylabel('CMF X'); legend; title('X Channel');

subplot(1,3,2);
plot(wl_peaks, cmf_band(:,2), 'go-', 'DisplayName','SPD-weighted Y'); hold on;
plot(wl_peaks, cmf_y,        'mx--', 'DisplayName','Peak-interp Y');
xlabel('Wavelength (nm)'); ylabel('CMF Y'); legend; title('Y Channel');

subplot(1,3,3);
plot(wl_peaks, cmf_band(:,3), 'co-', 'DisplayName','SPD-weighted Z'); hold on;
plot(wl_peaks, cmf_z,        'kx--', 'DisplayName','Peak-interp Z');
xlabel('Wavelength (nm)'); ylabel('CMF Z'); legend; title('Z Channel');
sgtitle('SPD-weighted CMFs vs Interpolated-at-Peak');

for i = 1:3
    subplot(1,3,i); grid on;
end
%% only CMFs vs CMFs x D50

figure;

subplot(3,1,1);
plot(wl_peaks, cmf_x, 'bo-', 'DisplayName', 'CMF_x'); hold on;
plot(wl_peaks, cmf_x .* d50_IP, 'ro--', 'DisplayName', 'CMF_x × D50');
xlabel('Wavelength (nm)'); ylabel('X');
title('X Channel: CMF vs. CMF × D50');
legend; grid on;

subplot(3,1,2);
plot(wl_peaks, cmf_y, 'go-', 'DisplayName', 'CMF_y'); hold on;
plot(wl_peaks, cmf_y .* d50_IP, 'mo--', 'DisplayName', 'CMF_y × D50');
xlabel('Wavelength (nm)'); ylabel('Y');
title('Y Channel: CMF vs. CMF × D50');
legend; grid on;

subplot(3,1,3);
plot(wl_peaks, cmf_z, 'co-', 'DisplayName', 'CMF_z'); hold on;
plot(wl_peaks, cmf_z .* d50_IP, 'ko--', 'DisplayName', 'CMF_z × D50');
xlabel('Wavelength (nm)'); ylabel('Z');
title('Z Channel: CMF vs. CMF × D50');
legend; grid on;

sgtitle('Peak-Interpolated CMFs vs. D50-Weighted CMFs');
%% normalised CMFs vs normalised CMFs x D50
k_CMFs     = sum(cmf_y);           
k_D50CMFs  = sum(d50_IP .* cmf_y); 

cmf_x_norm  = cmf_x  / k_CMFs;
cmf_y_norm  = cmf_y  / k_CMFs;
cmf_z_norm  = cmf_z  / k_CMFs;


cmf_xD50_norm = (cmf_x .* d50_IP) / k_D50CMFs;
cmf_yD50_norm = (cmf_y .* d50_IP) / k_D50CMFs;
cmf_zD50_norm = (cmf_z .* d50_IP) / k_D50CMFs;

figure;
subplot(3,1,1);
plot(wl_peaks, cmf_x_norm,      'bo-', 'DisplayName','CMF_x / k');
hold on;
plot(wl_peaks, cmf_xD50_norm,   'ro--', 'DisplayName','(CMF_x × D50) / k');
ylabel('X'); legend; grid on;
title('Normalized CMF_x: plain vs D50-weighted');

subplot(3,1,2);
plot(wl_peaks, cmf_y_norm,      'go-', 'DisplayName','CMF_y / k');
hold on;
plot(wl_peaks, cmf_yD50_norm,   'mo--', 'DisplayName','(CMF_y × D50) / k');
ylabel('Y'); legend; grid on;
title('Normalized CMF_y: plain vs D50-weighted');

subplot(3,1,3);
plot(wl_peaks, cmf_z_norm,      'co-', 'DisplayName','CMF_z / k');
hold on;
plot(wl_peaks, cmf_zD50_norm,   'ko--', 'DisplayName','(CMF_z × D50) / k');
xlabel('Wavelength (nm)');
ylabel('Z'); legend; grid on;
title('Normalized CMF_z: plain vs D50-weighted');

sgtitle('Normalized CMFs vs Normalized (CMFs × D50)');

%% normalised with a redish illumiannt

fake_illum = ones(size(wl_peaks));
fake_illum(wl_peaks < 550) = 0;  

k_CMFs      = sum(cmf_y);                    
k_FakeIllum = sum(fake_illum .* cmf_y);     

cmf_x_norm = cmf_x / k_CMFs;
cmf_y_norm = cmf_y / k_CMFs;
cmf_z_norm = cmf_z / k_CMFs;

cmf_xFake_norm = (cmf_x .* fake_illum) / k_FakeIllum;
cmf_yFake_norm = (cmf_y .* fake_illum) / k_FakeIllum;
cmf_zFake_norm = (cmf_z .* fake_illum) / k_FakeIllum;

figure;
subplot(3,1,1);
plot(wl_peaks, cmf_x_norm,     'bo-', 'DisplayName','CMF_x / k');
hold on;
plot(wl_peaks, cmf_xFake_norm, 'ro--', 'DisplayName','(CMF_x × FakeIllum) / k');
ylabel('X'); legend; grid on;
title('Normalized CMF_x: plain vs fake-illuminant');

subplot(3,1,2);
plot(wl_peaks, cmf_y_norm,     'go-', 'DisplayName','CMF_y / k');
hold on;
plot(wl_peaks, cmf_yFake_norm, 'mo--', 'DisplayName','(CMF_y × FakeIllum) / k');
ylabel('Y'); legend; grid on;
title('Normalized CMF_y: plain vs fake-illuminant');

subplot(3,1,3);
plot(wl_peaks, cmf_z_norm,     'co-', 'DisplayName','CMF_z / k');
hold on;
plot(wl_peaks, cmf_zFake_norm, 'ko--', 'DisplayName','(CMF_z × FakeIllum) / k');
xlabel('Wavelength (nm)');
ylabel('Z'); legend; grid on;
title('Normalized CMF_z: plain vs fake-illuminant');

sgtitle('Normalized CMFs vs Normalized (CMFs × Red Illuminant)');
