% process_film_scans.m
% Script to build, white-balance, and save multispectral film scans
% stored as TIFFs in "FlatFielded" subfolders, with interactive ROI cropping.

clear; close all;

%% User parameters
rootFolder = '/home/oem/eliza/data/film_scans/to_process';  
outputDir  = '/home/oem/eliza/data/film_scans/white_balanced';
led_path   = '../../../data/film/CVCL10bands.txt';
D50_path   = '../../../data/CIE_D50.txt';
cmf_path   = '../../../data/CIE2degCMFs_full.csv';

%% Read LED SPD and compute band center-wavelengths
LEDset   = readmatrix(led_path, 'Delimiter','\t');  
led_wl   = LEDset(:,1);            % fine sampling wavelengths
led_spds = LEDset(:,2:end);        % SPDs per band (length×B)
bandN    = size(led_spds,2);
wl = (sum(led_spds .* led_wl,1) ./ sum(led_spds,1))';  % center wavelengths

%% Precompute D50-weighted CMFs for rendering
D50    = importdata(D50_path);               % [λ, SPD]
CMFs   = importdata(cmf_path);               % [λ, x,y,z]
d50_spd = interp1(D50(:,1), D50(:,2), led_wl, 'linear','extrap');
cmf_xyz = interp1(CMFs(:,1), CMFs(:,2:4), led_wl, 'linear','extrap');

d50_band = zeros(bandN,1);
cmf_band = zeros(bandN,3);
for k = 1:bandN
    spd = led_spds(:,k);
    norm_spd = spd / sum(spd);
    d50_band(k)   = sum(d50_spd .* norm_spd);
    cmf_band(k,:) = sum(cmf_xyz .* norm_spd,1);
end
% normalization constant for XYZ→RGB
data_k_norm = sum(d50_band .* cmf_band(:,2));

%% Find all scan subfolders under rootFolder
D = dir(rootFolder);
isub = [D.isdir] & ~ismember({D.name},{'.','..'});
children = D(isub);

%% Process each scan
for i = 1:numel(children)
    childName = children(i).name;
    ffFolder  = fullfile(rootFolder, childName, 'FlatFielded');
    if ~exist(ffFolder,'dir')
        warning('Skipping %s: no FlatFielded folder.', childName);
        continue;
    end
    
    %% 1) Build the multispectral cube from TIFFs
    tiffs = dir(fullfile(ffFolder,'*.tif'));
    if isempty(tiffs)
        warning('No TIFFs in %s', ffFolder);
        continue;
    end
    img0 = imread(fullfile(ffFolder, tiffs(1).name));
    [H,W,~] = size(img0);
    cube_uint16 = zeros(H,W,numel(tiffs),'uint16');
    for k = 1:numel(tiffs)
        cube_uint16(:,:,k) = imread(fullfile(ffFolder, tiffs(k).name));
    end
    cube = double(cube_uint16) / double(intmax('uint16'));  % [0,1]
    
    %% 2) Select white patch via band 9
    band9_adj = imadjust(min(cube(:,:,9), 1));
    hW = figure; imshow(band9_adj); axis image off;
    title('Draw rectangle around white patch','FontSize',14);
    hR = drawrectangle('StripeColor','r'); wait(hR);
    pos = round(hR.Position);
    close(hW);
    x1 = max(1,pos(1)); y1 = max(1,pos(2));
    x2 = min(W, pos(1)+pos(3)-1);
    y2 = min(H, pos(2)+pos(4)-1);
    
    %% 3) Compute white_patch spectrum
    patchCube   = cube(y1:y2, x1:x2, :);
    white_patch = squeeze(mean(mean(patchCube,1),2));  % [B×1]
    
    %% 4) White-balance cube
    wb_cube_full = bsxfun(@times, cube, 1./reshape(white_patch,1,1,bandN));
    wb_cube_full(isnan(wb_cube_full)) = 0;
    wb_cube_full = min(max(wb_cube_full,0),1);
    
    %% 5) D50 rendering AFTER white balancing (sRGB)
    rad_wb      = wb_cube_full .* reshape(d50_band,1,1,bandN);
    lin_wb      = reshape(rad_wb,[],bandN) * cmf_band;   % [N×3]
    XYZ_wb_flat = lin_wb / data_k_norm;
    XYZimg2     = reshape(XYZ_wb_flat, H, W, 3);
    RGB_after   = xyz2rgb(XYZimg2, 'ColorSpace','srgb', 'WhitePoint','d50');
    RGB_after   = max(min(RGB_after,1),0);
    hFig = figure; imshow(RGB_after); title('sRGB under D50 - After White Balance','FontSize',14);
    pause(2); close(hFig);
    outAfter = fullfile(outputDir, [childName '_srgb_after.png']);
    imwrite(RGB_after, outAfter);

     %% 6) Crop final cube before saving
    hC2 = figure; imshow(wb_cube_full(:,:,9)); axis image off;
    title('Draw rectangle to crop final result','FontSize',14);
    hRect = drawrectangle('StripeColor','g'); wait(hRect);
    pos2 = round(hRect.Position);
    close(hC2);  
    x3 = max(1,pos2(1)); y3 = max(1,pos2(2));
    x4 = min(W, x3 + pos2(3)-1);
    y4 = min(H, y3 + pos2(4)-1);

    wb_cube = wb_cube_full(y3:y4, x3:x4, :);
    
    %% 7) Save cropped data
    if ~exist(outputDir,'dir'), mkdir(outputDir); end
    outMat = fullfile(outputDir, [childName '_balanced.mat']);
    save(outMat, 'wb_cube', 'wl', '-v7.3');
    fprintf('Saved cropped MAT: %s\n', outMat);
end
