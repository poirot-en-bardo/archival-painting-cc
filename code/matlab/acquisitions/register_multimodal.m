clear;
% --- Set file paths and parameters ---
D50_path     = '../../../data/CIE_D50.txt';
CMF_path     = '../../../data/CIE2degCMFs_1931.txt';
CMF_path_full= '../../../data/CIE2degCMFs_full.csv';
led_path     = '../../../data/film/CVCL10bands.txt';
matPath      = '/home/oem/eliza/data/film_scans/white_balanced/cactus_halogen_kodak_exp0_balanced.mat';
fixPath      = '/home/oem/eliza/data/processed/reflectance/registered/cactus_reflectance_full_before.hdr';
refBandMS    = 9;

% --- Load cubes and wavelengths ---
S       = load(matPath, 'wb_cube', 'wl');
movCube = flip(S.wb_cube, 1);
wl_mov  = S.wl;
hcFix   = hypercube(fixPath);
fixCube = gather(hcFix.DataCube);
wl_fix  = hcFix.Wavelength;

% --- Interactive crop of MSI cube ---
threshold = 1;
bandImg = double(movCube(:,:,refBandMS));
bandImg(bandImg > threshold) = threshold;
bandImg = (bandImg - min(bandImg(:))) / (max(bandImg(:)) - min(bandImg(:)));
bandImgAdj = imadjust(bandImg);

figure; imshow(bandImgAdj, []);
axis image off;
title('Draw ROI on MOVING (MSI) cube','FontSize',14);
hRec = drawrectangle('StripeColor','r');
wait(hRec);
pos = round(hRec.Position);
close;

movCrop = movCube( ...
    pos(2):(pos(2)+pos(4)-1), ...
    pos(1):(pos(1)+pos(3)-1), ...
    : );
fixCrop = fixCube; % Or crop HSI similarly if desired

% --- LED spectral info for MSI rendering ---
LEDset   = readmatrix(led_path, 'Delimiter','\t');  
led_wl   = LEDset(:,1);             
led_spds = LEDset(:,2:end);         
bandN    = size(led_spds,2);

% --- Render HSI (fixed) cube to D50 sRGB ---
ill = importdata(D50_path);          
CMFs = importdata(CMF_path); 
illIP = interp1(ill(:,1),     ill(:,2), wl_fix, 'spline');        
CMFx  = interp1(CMFs(:,1), CMFs(:,2), wl_fix, 'spline');         
CMFy  = interp1(CMFs(:,1), CMFs(:,3), wl_fix, 'spline');
CMFz  = interp1(CMFs(:,1), CMFs(:,4), wl_fix, 'spline');
CMFsIP = [CMFx, CMFy, CMFz];
S_fix = CMFsIP .* illIP; % [B_fix x 3]
k_fix = sum(S_fix(:,2));
linCube_fix = reshape(fixCrop, [], length(wl_fix));
XYZ_fix = (linCube_fix * S_fix) ./ k_fix;
XYZimg_fix = reshape(XYZ_fix, size(fixCrop,1), size(fixCrop,2), 3);
RGB_fix = xyz2rgb(XYZimg_fix, 'ColorSpace','srgb', 'WhitePoint','d50');
RGB_fix = max(min(RGB_fix,1),0);

% --- Render MSI (moving) cube to D50 sRGB (band-averaged) ---
D50      = importdata(D50_path);
CMFs_f   = importdata(CMF_path_full);
d50_spd  = interp1(D50(:,1), D50(:,2), led_wl, 'linear','extrap');
cmf_xyz  = interp1(CMFs_f(:,1), CMFs_f(:,2:4), led_wl, 'linear','extrap');
d50_band = zeros(bandN,1);
cmf_band = zeros(bandN,3);
for k = 1:bandN
    spd = led_spds(:,k);
    norm_spd = spd / sum(spd);
    d50_band(k)   = sum(d50_spd .* norm_spd);
    cmf_band(k,:) = sum(cmf_xyz .* norm_spd,1);
end
data_k_norm = sum(d50_band .* cmf_band(:,2));
rad      = movCrop .* reshape(d50_band,1,1,bandN);
linCube  = reshape(rad,[],bandN) * cmf_band;
XYZ_mov  = linCube / data_k_norm;
XYZimg_mov = reshape(XYZ_mov, size(movCrop,1), size(movCrop,2), 3);
RGB_mov = xyz2rgb(XYZimg_mov, 'ColorSpace','srgb', 'WhitePoint','d50');
RGB_mov = max(min(RGB_mov,1),0);
%%
% --- Register MSI RGB to HSI RGB (use grayscale images for registration) ---
gray_mov = imadjust(rgb2gray(RGB_mov));
gray_fix = imadjust(rgb2gray(RGB_fix));

[optimizer, metric] = imregconfig('multimodal');
optimizer.MaximumIterations = 500;
optimizer.InitialRadius = 1e-3;

tformObj = imregtform(gray_mov, gray_fix, 'affine', optimizer, metric);

outputView = imref2d(size(gray_fix));
reg_mov_rgb = imwarp(RGB_mov, tformObj, 'OutputView', outputView);

% --- Show overlay AFTER registration only ---
figure; imshowpair(RGB_fix, reg_mov_rgb, 'montage');
title('Overlay AFTER Registration (blend)','FontSize',14);
%%
figure; imshowpair(RGB_fix, reg_mov_rgb);
title('D50 RGBs AFTER Registration (Left: Fixed HSI, Right: Reg. MSI)','FontSize',14);
%%

% --- Apply transform to every band of cropped, flipped MSI cube ---
[rowsFix, colsFix] = size(RGB_fix(:,:,1));
registeredCube = zeros(rowsFix, colsFix, bandN, 'like', movCrop);
for b = 1:bandN
    band = movCrop(:,:,b);
    registeredCube(:,:,b) = imwarp(band, tformObj, 'OutputView', outputView);
end

disp('Done! Overlay and registered MSI cube are available.');
