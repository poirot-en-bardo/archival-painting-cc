clear; close all;
input_folder = '/home/oem/eliza/mac-shared/registered';  % Where old .mat files are
output_folder = fullfile(input_folder, 'reshaped');

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

mat_files = dir(fullfile(input_folder, '*.mat'));

for i = 1:length(mat_files)
    file_path = fullfile(input_folder, mat_files(i).name);
    fprintf('Processing: %s\n', mat_files(i).name);
    
    data = load(file_path);

    % Check that needed fields exist
    if ~all(isfield(data, {'mov_XYZ', 'mov_Lab', 'mov_RGB'}))
        warning('Missing expected fields in %s. Skipping.', mat_files(i).name);
        continue;
    end

    % Infer dimensions from number of pixels
    nPixels = size(data.mov_XYZ, 1);
    W = sqrt(nPixels); H = nPixels / W;
    if abs(H - round(H)) > 0 || abs(W - round(W)) > 0
        warning('Cannot infer square shape in %s. Skipping.', mat_files(i).name);
        continue;
    end
    H = round(H); W = round(W);

    % Reshape
    XYZ_img = reshape(data.mov_XYZ, H, W, 3);
    Lab_img = reshape(data.mov_Lab, H, W, 3);
    RGB_img = reshape(data.mov_RGB, H, W, 3);

    % Save reshaped copy to output folder
    save_path = fullfile(output_folder, mat_files(i).name);
    save(save_path, 'XYZ_img', 'Lab_img', 'RGB_img');
    fprintf('Saved reshaped file to: %s\n', save_path);
end
