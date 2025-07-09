% Set folders
src_folder = '/home/oem/eliza/data/xyz_lab_rgb/colorchecker';
dst_folder = '/home/oem/eliza/data/xyz_lab_rgb/colorchecker/recomputed';
if ~exist(dst_folder, 'dir')
    mkdir(dst_folder);
end

% Get list of .mat files
files = dir(fullfile(src_folder, '*.mat'));

% Loop through each .mat file
for k = 1:length(files)
    fprintf('Processing %s...\n', files(k).name);
    
    % Load data
    src_path = fullfile(src_folder, files(k).name);
    data = load(src_path);
    
    % Extract and normalize patch XYZ
    XYZ = data.patchXYZ;  % [24 x 3]
    XYZ_normalized = XYZ ./ 100;

    % Recompute Lab
    data.patchLab = xyz2lab_custom(XYZ);

    % Recompute RGB and RGB_lin
    data.patchRGB_lin = xyz2prophoto(XYZ_normalized, false);
    data.patchRGB     = xyz2prophoto(XYZ_normalized, true);

    % Optional: Visualize RGB swatches
    rgb_vis = reshape(data.patchRGB, [1, 24, 3]);
    figure('Name', ['Patch RGB: ' files(k).name], 'Color', 'w');
    imshow(rgb_vis, 'InitialMagnification', 'fit');
    title(['RGB Patch Preview: ' files(k).name], 'Interpreter', 'none', 'FontSize', 14);
    drawnow; pause(0.3);

    % Save updated file
    dst_path = fullfile(dst_folder, files(k).name);
    save(dst_path, '-struct', 'data', '-v7.3');
end
