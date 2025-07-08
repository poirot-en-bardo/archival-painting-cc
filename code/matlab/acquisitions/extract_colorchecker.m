inputFile = '/home/oem/eliza/mac-shared/to_register/full/cactus_after.mat';  
data = load(inputFile);

% Check required fields
if ~isfield(data, 'RGB_img') || ~isfield(data, 'XYZ_img') || ~isfield(data, 'Lab_img')
    error('MAT file must contain RGB_img, XYZ_img, and Lab_img.');
end

disp_img = im2uint8(data.RGB_img);  % Increase brightness by 50%
disp_img = min(disp_img, 255); 

% Display RGB image
figure('Name','Select ColorChecker Patches');
imshow(disp_img);
title('Click and drag to draw rectangles around 24 patches (top-left to bottom-right)');

% Initialize storage
numPatches = 24;
patchRGB = zeros(numPatches, 3);
patchXYZ = zeros(numPatches, 3);
patchLab = zeros(numPatches, 3);
positions = zeros(numPatches, 4);

for i = 1:numPatches
    h = drawrectangle('Label', sprintf('%d', i), 'Color', 'r');
    wait(h);
    pos = round(h.Position);
    positions(i,:) = pos;
    
    % Crop each region (with bounds check)
    x1 = max(1, pos(1)); y1 = max(1, pos(2));
    x2 = min(size(data.RGB_img,2), x1 + pos(3) - 1);
    y2 = min(size(data.RGB_img,1), y1 + pos(4) - 1);
    
    rgb_crop = data.RGB_img(y1:y2, x1:x2, :);
    xyz_crop = data.XYZ_img(y1:y2, x1:x2, :);
    lab_crop = data.Lab_img(y1:y2, x1:x2, :);

    % Store mean patch values
    patchRGB(i,:) = squeeze(mean(reshape(rgb_crop, [], 3), 1));
    patchXYZ(i,:) = squeeze(mean(reshape(xyz_crop, [], 3), 1));
    patchLab(i,:) = squeeze(mean(reshape(lab_crop, [], 3), 1));
end

% Save results in same folder with _colorchecker suffix
outputDir = '/home/oem/eliza/mac-shared/colorchecker';
[path, name, ~] = fileparts(inputFile);
save_name = fullfile(outputDir, [name '_colorchecker.mat']);
save(save_name, 'patchRGB', 'patchXYZ', 'patchLab', 'positions');

fprintf('Saved color checker patch data to: %s\n', save_name);
