

clear; close all; clc;

img_path1 = '/Volumes/School/Thesis/data/captures/registered/prophoto/yoda1_ref_hsi_after.png';
img_path2 = '/Volumes/School/Thesis/data/captures/registered/prophoto/yoda1_reg_hsi_before.png';
img_path1= '/Volumes/School/Thesis/data/captures/registered/prophoto/yoda2_ref_hsi_kodak_halogen_after.png';
img_path2 = '/Volumes/School/Thesis/data/captures/registered/prophoto/yoda2_reg_kodak_halogen_after.png';
saveDir  = 'results/wavelet_diff';
waveletName = 'db2';
numLevels   = 3;

%% 1) Read and preprocess images
I1_rgb = imread(img_path1);
I2_rgb = imread(img_path2);
I1 = im2double(rgb2gray(I1_rgb));
I2 = im2double(rgb2gray(I2_rgb));

%% 2) Create save directory if needed
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% 3) Perform multi-level wavelet decomposition
[C1,S1] = wavedec2(I1, numLevels, waveletName);
[C2,S2] = wavedec2(I2, numLevels, waveletName);

%% 4) Compute approximation (low-frequency) band difference
A1 = appcoef2(C1, S1, waveletName, numLevels);
A2 = appcoef2(C2, S2, waveletName, numLevels);
diffApprox = mat2gray(abs(A1 - A2));
% imwrite(diffApprox, fullfile(saveDir, 'approx_diff.png'));

disp('Saved approximation band difference: approx_diff.png');

%% 5) Compute and save detail (high-frequency) band differences
for lvl = 1:numLevels
    [H1, V1, D1] = detcoef2('all', C1, S1, lvl);
    [H2, V2, D2] = detcoef2('all', C2, S2, lvl);
    diffH = mat2gray(abs(H1 - H2));
    diffV = mat2gray(abs(V1 - V2));
    diffD = mat2gray(abs(D1 - D2));
    % imwrite(diffH, fullfile(saveDir, sprintf('level%d_horizontal_diff.png', lvl)));
    % imwrite(diffV, fullfile(saveDir, sprintf('level%d_vertical_diff.png',   lvl)));
    % imwrite(diffD, fullfile(saveDir, sprintf('level%d_diagonal_diff.png',   lvl)));
    fprintf('Saved detail diffs for level %d: horizontal, vertical, diagonal\n', lvl);
end

%% 6) (Optional) Display summary of sub-bands
figure('Name','Wavelet Difference Summary','NumberTitle','off','Position',[100 100 1200 600]);
subplot(1, numLevels+1, 1), imshow(I1), title('Image 1');
subplot(1, numLevels+1, 2), imshow(I2), title('Image 2');
subplot(1, numLevels+1, 3), imshow(diffApprox), title('Approx Diff');
% for lvl = 1:numLevels
%     idxRow1 = 1;
%     idxRow2 = 2;
%     idxCol = lvl + 2;
%     % Show detail from image1 for visual reference (optional)
%     subplot(2, numLevels+1, (idxRow1-1)*(numLevels+1) + idxCol);
%     imshow(mat2gray(detcoef2('h',C1,S1,lvl)));
%     title(sprintf('L%d Orig H', lvl));
%     subplot(2, numLevels+1, (idxRow2-1)*(numLevels+1) + idxCol);
%     imshow(mat2gray(detcoef2('h',C2,S2,lvl)));
%     title(sprintf('L%d Reg H', lvl));
% end


disp('Wavelet-based difference analysis complete.');
%%
% Option 1: Manual threshold
thresh = 0.1; % Try values between 0.05 and 0.2 depending on sensitivity
mask = diffApprox > thresh;

% Option 2: Otsu's method (automatic threshold)
level = graythresh(diffApprox);
mask_otsu = diffApprox > level;

% Save or show the mask
figure; imshow(mask); title('Change Mask (Manual Threshold)');
figure; imshow(mask_otsu); title('Change Mask (Otsu Threshold)');