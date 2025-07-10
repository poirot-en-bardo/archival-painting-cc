% --- Parameters ---
binSize = 2;
inputPath = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat';  % Replace with actual path

data = load(inputPath);
fields = fieldnames(data);

% --- Extract and bin spec_mask ---
spec_mask_orig = data.spec_mask;

[H, W] = size(spec_mask_orig);
h = floor(H / binSize);
w = floor(W / binSize);

% Crop and reshape for binning
spec_mask_crop = spec_mask_orig(1:h*binSize, 1:w*binSize);
spec_mask_crop = reshape(spec_mask_crop, binSize, h, binSize, w);
spec_mask_crop = permute(spec_mask_crop, [2 4 1 3]);  % [h w binSize binSize]

% Average and threshold
spec_mask_binned = mean(mean(spec_mask_crop, 3), 4) > 0.5;

% --- Replace original field with binned version ---
data.spec_mask = spec_mask_binned;

% --- Save to new file in same folder ---
[inputFolder, baseName, ~] = fileparts(inputPath);
outputPath = fullfile(inputFolder, [baseName '_fixed.mat']);
save(outputPath, '-struct', 'data');

fprintf('Saved fixed file to: %s\n', outputPath);