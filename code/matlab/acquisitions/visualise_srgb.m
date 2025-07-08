clear; close all;

% --- Input and output folders ---
inputFolder  = '/home/oem/eliza/data/xyz_lab_rgb/colorchecker';    
outputFolder = '/home/oem/eliza/data/xyz_lab_rgb/srgb_cc'; 

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% --- List all .mat files ---
files = dir(fullfile(inputFolder, '*.mat'));

for i = 1:length(files)
    matPath = fullfile(files(i).folder, files(i).name);
    fprintf('Processing: %s\n', files(i).name);
    
    % --- Load .mat file ---
    S = load(matPath);
    % if ~isfield(S, 'XYZ_img')
    %     warning('Skipping %s: No field named "XYZ_img".', files(i).name);
    %     continue;
    % end

    % XYZ = S.XYZ_img;
    % 
    % % --- Normalize to 0–1 range for sRGB conversion ---
    % XYZ_norm = XYZ ./ 100;  % assumes XYZ in [0–100]
    % 
    % % --- Convert to sRGB (D50 reference) ---
    % sRGB = xyz2rgb(XYZ_norm, 'ColorSpace', 'srgb', 'WhitePoint', 'd50');
    % sRGB = min(max(sRGB, 0), 1);  % Clamp to [0,1]
    % 
    % % --- Save image ---
    % [~, baseName, ~] = fileparts(files(i).name);
    % outPath = fullfile(outputFolder, [baseName '_srgb.png']);
    % imwrite(sRGB, outPath);
    % 
    % fprintf('Saved: %s\n', outPath);
    %% --- Visualize the 6×4 grid of selected patch RGBs
    patchXYZ = S.patchXYZ;
    patchRGB_srgb = xyz2rgb(patchXYZ ./100, 'ColorSpace', 'srgb', 'WhitePoint', 'd50');
    patchRGB_srgb = min(max(patchRGB_srgb, 0), 1);  % Clamp to [0,1]
    % patchRGB_srgb = patchRGB;
    % Parameters
    rows = 4;
    cols = 6;
    patchSize = 70;  % Size of each square patch (pixels)
    numPatches = 24;
    % Create grid image
    gridImg = zeros(patchSize * rows, patchSize * cols, 3);
    
    for j = 1:numPatches
        r = floor((j-1)/cols);  % row index
        c = mod((j-1), cols);   % col index
    
        color = reshape(patchRGB_srgb(j,:), 1, 1, 3);
        block = repmat(color, patchSize, patchSize);
    
        rowStart = r * patchSize + 1;
        colStart = c * patchSize + 1;
        gridImg(rowStart:rowStart+patchSize-1, colStart:colStart+patchSize-1, :) = block;
    end
    
    sRGB = min(max(gridImg, 0), 1);  % Just in case

    [~, baseName, ~] = fileparts(files(i).name);
    outPath = fullfile(outputFolder, [baseName '_srgb.png']);
    imwrite(sRGB, outPath);

end

fprintf('All files processed.\n');
