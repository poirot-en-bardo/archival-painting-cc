
% Simplified pipeline: compute XYZ/Lab/RGB only for FlatFielded folders
clear; close all;
parent_folder = "/home/oem/eliza/data/to_register/cactus";
outFolder = '/home/oem/eliza/mac-shared/registered';
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% Load CMFs and D50
cmf_path = '../../../data/CIE2degCMFs_full.csv';
D50_path = '../../../data/CIE_D50.txt';
led_path = '../../../data/film/CVCL10bands.txt';

fullCMFs = importdata(cmf_path);
D50_SPD = importdata(D50_path);
LEDset = readmatrix(led_path, 'Delimiter','\t');
wl_led = LEDset(:,1);
bandN = size(LEDset, 2) - 1;

% Setup folder list
folders = dir(parent_folder);
folders = folders([folders.isdir]);
folders = folders(~ismember({folders.name},{'.','..'}));

scale = 0.4;
cube_idx = 1;

for i = 1:length(folders)
    mov_folder = fullfile(parent_folder, folders(i).name);
    flatfielded = fullfile(mov_folder, 'FlatFielded');
    if ~exist(flatfielded, 'dir')
        fprintf('Skipping %s: no FlatFielded folder.\n', folders(i).name);
        continue;
    end

    fprintf('Processing %s...\n', folders(i).name);

    tiff_files = dir(fullfile(flatfielded, '*.tif'));
    if isempty(tiff_files)
        warning('No valid TIFF image data in %s', flatfielded);
        continue;
    end

    sample = imread(fullfile(flatfielded, tiff_files(1).name));
    [H, W] = size(sample);
    cube_mov = zeros(H, W, numel(tiff_files), class(sample));
    for k = 1:numel(tiff_files)
        cube_mov(:,:,k) = imread(fullfile(flatfielded, tiff_files(k).name));
    end

    cube_mov = imresize(cube_mov, scale);
    cube_mov = flip(cube_mov, 2);  % Flip horizontally
    cube_mov = flip(cube_mov, 1);  % Flip vertically

    slice = double(cube_mov(:,:,round(size(cube_mov,3)/2)));
    slice = slice - min(slice(:));
    slice = slice / max(slice(:));
    slice = slice .^ 0.2;  % high gamma for visibility
    figure; imagesc(slice); axis image; colormap gray; title('Crop TIFF image');
    roi = drawrectangle('InteractionsAllowed','all'); wait(roi); 
    r = round(roi.Position);
    xa = max(1, r(1)); ya = max(1, r(2));
    xb = min(size(cube_mov,2), xa + r(3) - 1); yb = min(size(cube_mov,1), ya + r(4) - 1);
    cube_mov = cube_mov(ya:yb, xa:xb, :);

    d50_spd = interp1(D50_SPD(:,1), D50_SPD(:,2), wl_led, 'linear', 'extrap');
    cmf_xyz = interp1(fullCMFs(:,1), fullCMFs(:,2:4), wl_led, 'linear', 'extrap');

    d50_band = zeros(bandN, 1);
    cmf_band = zeros(bandN, 3);
    for k = 1:bandN
        led_spd = LEDset(:, k+1);
        norm_led = led_spd / sum(led_spd);
        d50_band(k) = sum(d50_spd .* norm_led);
        cmf_band(k, :) = sum(cmf_xyz .* norm_led, 1);
    end

    cube_dbl = double(cube_mov) / double(intmax('uint16'));
    refl = reshape(cube_dbl, [], bandN);
    mov_XYZ = ref2xyz(d50_band(:), cmf_band, refl);


    mov_Lab = xyz2lab(mov_XYZ);
    mov_RGB = xyz2prophoto(mov_XYZ ./ 100, true);

    H = size(cube_mov, 1);
    W = size(cube_mov, 2);
    XYZ_img = reshape(mov_XYZ, H, W, 3);
    Lab_img = reshape(mov_Lab, H, W, 3);
    RGB_img = reshape(mov_RGB, H, W, 3);

    figure; imagesc(mat2gray(XYZ_img)); axis image; title('XYZ Image Preview');
    figure; imagesc(mat2gray(RGB_img)); axis image; title('RGB Image Preview');

    save_name = sprintf('%s_data_%d', folders(i).name, cube_idx);
    save(fullfile(outFolder, [save_name '.mat']), 'XYZ_img', 'Lab_img', 'RGB_img');
    cube_idx = cube_idx + 1;
end
