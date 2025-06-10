% Hyperspectral XYZ/Lab/RGB computation pipeline
clear; close all;
parent_folder = "/home/oem/eliza/data/to_register/yoda";
outFolder = '/home/oem/eliza/mac-shared/registered_hyper';
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% Load CMFs and D50
cmf_path = '../../../data/CIE2degCMFs_full.csv';
D50_path = '../../../data/CIE_D50.txt';

fullCMFs = importdata(cmf_path);
D50_SPD = importdata(D50_path);

% Setup folder list
folders = dir(fullfile(parent_folder, '**/')); % recursive
folders = folders([folders.isdir]);
folders = folders(~ismember({folders.name},{'.','..'}));

scale = 0.4;
cube_idx = 1;

for i = 1:length(folders)
    mov_folder = fullfile(folders(i).folder, folders(i).name);
    fprintf('Looking in: %s\n', mov_folder);

    hdr_file = dir(fullfile(mov_folder, '*.hdr'));
    img_file = dir(fullfile(mov_folder, '*.img'));

    if isempty(hdr_file) || isempty(img_file)
        fprintf('Skipping %s: no hyperspectral cube found.\n', folders(i).name);
        continue;
    end

    fprintf('Processing %s...\n', folders(i).name);
    hcube = hypercube(fullfile(mov_folder, hdr_file(1).name));
    cube_mov = hcube.DataCube;
    wl_mov = hcube.Wavelength;
    valid_idx = find(wl_mov >= 380 & wl_mov <= 780);
    cube_mov = cube_mov(:,:,valid_idx);
    wl_mov = wl_mov(valid_idx);

    cube_mov = flip(cube_mov, 2);  % Flip horizontally
    cube_mov = imresize(cube_mov, scale);

    slice = double(cube_mov(:,:,round(size(cube_mov,3)/2)));
    slice = slice - min(slice(:));
    slice = slice / max(slice(:));
    slice = slice .^ 0.2;  % high gamma for visibility
    figure; imagesc(slice); axis image; colormap gray; title('Crop hyperspectral image');
    roi = drawrectangle('InteractionsAllowed','all'); wait(roi); 
    r = round(roi.Position);
    xa = max(1, r(1)); ya = max(1, r(2));
    xb = min(size(cube_mov,2), xa + r(3) - 1); yb = min(size(cube_mov,1), ya + r(4) - 1);
    cube_mov = cube_mov(ya:yb, xa:xb, :);

    refl = reshape(cube_mov, [], size(cube_mov,3));
    cmf_interp = interp1(fullCMFs(:,1), fullCMFs(:,2:4), wl_mov, 'linear', 'extrap');
    ill_interp = interp1(D50_SPD(:,1), D50_SPD(:,2), wl_mov, 'linear', 'extrap');
    mov_XYZ = ref2xyz(ill_interp(:), cmf_interp, refl);

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
