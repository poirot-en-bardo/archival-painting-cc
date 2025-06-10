% Main pipeline for hyperspectral/multispectral registration and colorimetric conversion
clear; close all;
parent_folder = "/home/oem/eliza/data/to_register/cactus";
reference_img_folder_name = "after";
outFolder = '/home/oem/eliza/mac-shared/registered';

% Load CMFs and D50
cmf = importdata('../../../data/CIE2degCMFs_1931.txt');
ill = importdata('../../../data/CIE_D50.txt');

% Setup folder list
folders = dir(parent_folder);
folders = folders([folders.isdir]);
folders = folders(~ismember({folders.name},{'.','..'}));

% Find reference "after" cube
ref_folder = fullfile(parent_folder, reference_img_folder_name);
ref_hdr = dir(fullfile(ref_folder, '*.hdr'));
ref_img = dir(fullfile(ref_folder, '*.img'));

if isempty(ref_hdr) || isempty(ref_img)
    error('Reference folder must contain .hdr and .img hypercube files.');
end

ref_hcube = hypercube(fullfile(ref_folder, ref_hdr(1).name));
cube_ref = ref_hcube.DataCube;
wl_ref = ref_hcube.Wavelength;
valid_idx = find(wl_ref >= 380 & wl_ref <= 780);
cube_ref = cube_ref(:,:,valid_idx);

% Flip horizontally
cube_ref = flip(cube_ref, 2);
% Subsample reference cube
scale = 0.4;
[nr, nc, nb] = size(cube_ref);
cube_ref_resized = zeros(round(nr*scale), round(nc*scale), nb);
cube_ref_resized = imresize(cube_ref, scale);
cube_ref = cube_ref_resized;

wl_ref = wl_ref(valid_idx);

ref_band = 80;

% User crop on reference
fig_ref = figure; imagesc(cube_ref(:,:,ref_band)); axis image; colormap gray; title('Crop reference image');
roi_ref = drawrectangle('InteractionsAllowed','all');
uiwait(fig_ref);  % Wait for manual confirmation by closing the figure
r_ref = round(get(roi_ref, 'Position'));
x_ra = max(1, r_ref(1)); y_ra = max(1, r_ref(2));
x_rb = min(size(cube_ref,2), x_ra + r_ref(3) - 1); y_rb = min(size(cube_ref,1), y_ra + r_ref(4) - 1);
cube_ref = cube_ref(y_ra:y_rb, x_ra:x_rb, :);

cube_idx = 1;
%%
for i = 1:length(folders)
    if strcmp(folders(i).name, reference_img_folder_name)
        continue;
    end

    mov_folder = fullfile(parent_folder, folders(i).name);
    flatfielded = fullfile(mov_folder, 'FlatFielded');
    if ~exist(flatfielded, 'dir')
        continue;
    end

    fprintf('Processing %s...\n', folders(i).name);

    hdr_file = dir(fullfile(flatfielded, '*.hdr'));
    img_file = dir(fullfile(flatfielded, '*.img'));

    if ~isempty(hdr_file) && ~isempty(img_file)
        hcube = hypercube(fullfile(flatfielded, hdr_file(1).name));
        cube_mov = hcube.DataCube;
        wl_mov = hcube.Wavelength;
        valid_idx = find(wl_mov >= 380 & wl_mov <= 780);
        cube_mov = cube_mov(:,:,valid_idx);

        % Flip horizontally
        cube_mov = flip(cube_mov, 2);
% Subsample moving hypercube
[nr, nc, nb] = size(cube_mov);
cube_mov_resized = zeros(round(nr*scale), round(nc*scale), nb);
cube_mov_resized = imresize(cube_mov, scale);
cube_mov = cube_mov_resized;

        wl_mov = wl_mov(valid_idx);
        mov_band = 80;
    else
        tiff_files = dir(fullfile(flatfielded, '*.tif'));
        if isempty(tiff_files)
            warning('No valid image data in %s', flatfielded);
            continue;
        end
        sample = imread(fullfile(flatfielded, tiff_files(1).name));
        [H, W] = size(sample);
        cube_mov = zeros(H, W, numel(tiff_files), class(sample));
% Will be subsampled after loading
        for k = 1:numel(tiff_files)
            cube_mov(:,:,k) = imread(fullfile(flatfielded, tiff_files(k).name));
        end
        mov_band = 8;
% Subsample TIFF-based multispectral cube
[nr, nc, nb] = size(cube_mov);
cube_mov_resized = zeros(round(nr*scale), round(nc*scale), nb);
for b = 1:nb
    cube_mov_resized(:,:,b) = imresize(cube_mov(:,:,b), scale);
end
cube_mov = cube_mov_resized;
    end

    % User crop on moving image
fig_mov = figure; imagesc(cube_mov(:,:,mov_band)); axis image; colormap gray; title('Crop moving image');
roi_mov = drawrectangle('InteractionsAllowed','all');
uiwait(fig_mov);
    r_mov = round(get(roi_mov, 'Position'));
    x_ma = max(1, r_mov(1)); y_ma = max(1, r_mov(2));
    x_mb = min(size(cube_mov,2), x_ma + r_mov(3) - 1); y_mb = min(size(cube_mov,1), y_ma + r_mov(4) - 1);
    cube_mov = cube_mov(y_ma:y_mb, x_ma:x_mb, :);

    I1 = mat2gray(cube_ref(:,:,ref_band));
    I2 = mat2gray(cube_mov(:,:,mov_band));

    % Try registration with user loop
    maxRatio = 0.6;
    registered = false;
    while ~registered
        points1 = detectSURFFeatures(I1);
        points2 = detectSURFFeatures(I2);
        [features1, v1] = extractFeatures(I1, points1);
        [features2, v2] = extractFeatures(I2, points2);
        pairs = matchFeatures(features1, features2, 'MaxRatio', maxRatio, 'Unique', true);

        mp1 = v1(pairs(:,1));
        mp2 = v2(pairs(:,2));
        [tform, inlierIdx] = estimateGeometricTransform2D(mp2, mp1, 'affine');

        figure; showMatchedFeatures(I1, I2, mp1, mp2, 'montage');
        title(sprintf('Matched points (MaxRatio=%.2f)', maxRatio));

        figure; imshowpair(I1, imwarp(I2, tform, 'OutputView', imref2d(size(I1))));
        title('Registered overlay');

        user_input = input('Is the registration good? (y/n): ', 's');
        if strcmpi(user_input, 'y')
            registered = true;
        else
            maxRatio = maxRatio * 0.9;
            fprintf('Reducing MaxRatio to %.2f\n', maxRatio);
        end
    end

    Rfixed = imref2d(size(I1));
    reg_mov = zeros(size(cube_ref));
    for b = 1:size(cube_ref,3)
        reg_mov(:,:,b) = imwarp(cube_mov(:,:,b), tform, 'OutputView', Rfixed);
    end

    valid = any(cube_ref > 0.01, 3) & any(reg_mov > 0.01, 3);
    [rows, cols] = find(valid);
    cube_ref_crop = cube_ref(min(rows):max(rows), min(cols):max(cols), :);
    reg_mov_crop = reg_mov(min(rows):max(rows), min(cols):max(cols), :);

    % Spectral to color conversion
    cmf_interp = interp1(cmf(:,1), cmf(:,2:4), wl_ref, 'linear', 'extrap');
    ill_interp = interp1(ill(:,1), ill(:,2), wl_ref, 'linear', 'extrap');

    ref_XYZ = ref2xyz(ill_interp(:), cmf_interp, reshape(cube_ref_crop, [], size(cube_ref_crop,3)));
    mov_XYZ = ref2xyz(ill_interp(:), cmf_interp, reshape(reg_mov_crop, [], size(reg_mov_crop,3)));

    ref_Lab = xyz2lab(ref_XYZ);
    mov_Lab = xyz2lab(mov_XYZ);

    ref_RGB = xyz2prophoto(ref_XYZ / 100, true);
    mov_RGB = xyz2prophoto(mov_XYZ / 100, true);

    % Save results
    save_name = sprintf('%s_reg', folders(i).name);
    save(fullfile(outFolder, [save_name '_ref_' num2str(cube_idx) '.mat']), 'ref_XYZ', 'ref_Lab', 'ref_RGB');
    save(fullfile(outFolder, [save_name '_mov_' num2str(cube_idx) '.mat']), 'mov_XYZ', 'mov_Lab', 'mov_RGB');
    cube_idx = cube_idx + 1;
end
