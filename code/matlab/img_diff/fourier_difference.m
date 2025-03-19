
img_folder = '../../../data/scream/mod';
% mod_img_path = '../../../data/scream/mod/scream_low_exp.png';
mod_img_path = '../../../results/illum_corrected1.png';
img_files = dir(fullfile(img_folder, '*.png'));  % Adjust for your file type
reference_path = '../../../data/scream/scream_rgb.png';
save_folder = '../../../results/fourier';

% For one image
process_fourier(reference_path, mod_img_path, save_folder);

% For a folder of images
% for i = 1:length(img_files)
%     img_path1 = fullfile(img_folder, img_files(i).name);
%     process_fourier(img_path1, reference_path, save_folder);
% end

function process_fourier(img_path1, img_path2, save_folder)
    % Ensure save folder exists
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end

    % Max value for uint16 normalization
    roof = double(intmax('uint16'));

    % Read images
    img1 = imread(img_path1);
    img2 = imread(img_path2);

    % Convert images to grayscale
    img1_grey = rgb2gray(img1);
    img2_grey = rgb2gray(img2);

    % Normalize to [0, 1]
    img1_grey = double(img1_grey) / roof;
    img2_grey = double(img2_grey) / roof;

    % Adjust img2 to match the mean intensity of img1 (exposure normalization)
    mean1 = mean(img1_grey(:));
    mean2 = mean(img2_grey(:));
    img2_grey = img2_grey * (mean1 / mean2);

    % Histogram matching for exposure normalization
    % img2_grey = imhistmatch(img2_grey, img1_grey);

    % Fourier Transform
    F1 = fft2(img1_grey);
    F2 = fft2(img2_grey);

    % Shift zero frequency to center
    F1_shifted = fftshift(F1);
    F2_shifted = fftshift(F2);

    % Extract frequencies
    [rows, cols] = size(F1_shifted);
    low_freq_radius = min(rows, cols) / 500;
    [X, Y] = meshgrid(1:cols, 1:rows);
    centerX = floor(cols / 2);
    centerY = floor(rows / 2);
    distance = sqrt((X - centerX).^2 + (Y - centerY).^2);
    low_freq_mask = distance <= low_freq_radius;
    high_freq_mask = distance > low_freq_radius;

    % Apply masks
    low_freq1 = F1_shifted .* low_freq_mask;
    high_freq1 = F1_shifted .* high_freq_mask;
    low_freq2 = F2_shifted .* low_freq_mask;
    high_freq2 = F2_shifted .* high_freq_mask;

    % Inverse FFT
    low_freq1_spatial = abs(ifft2(ifftshift(low_freq1)));
    high_freq1_spatial = abs(ifft2(ifftshift(high_freq1)));
    low_freq2_spatial = abs(ifft2(ifftshift(low_freq2)));
    high_freq2_spatial = abs(ifft2(ifftshift(high_freq2)));

    % Differences
    low_freq_diff = abs(low_freq1_spatial - low_freq2_spatial);
    high_freq_diff = abs(high_freq1_spatial - high_freq2_spatial);

    % Normalize for display
    low_freq_diff_norm = mat2gray(low_freq_diff);
    high_freq_diff_norm = mat2gray(high_freq_diff);

    % Save figures
    [~, img_name, ~] = fileparts(img_path2);
    
    figure; imshow(low_freq_diff_norm, []); colorbar; title(['Low Frequency Difference - ',img_name], Interpreter="none");
    saveas(gcf, fullfile(save_folder, [img_name '_low_freq_diff.png']));
    % close;

    figure; imshow(high_freq_diff_norm, []); colorbar;  title(['High Frequency Difference - ',img_name], Interpreter="none");
    saveas(gcf, fullfile(save_folder, [img_name '_high_freq_diff.png']));
    % close;
end
