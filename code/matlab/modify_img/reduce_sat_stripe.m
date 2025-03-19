% MATLAB script to locally reduce saturation and contrast in a selected region of the image

% Read the PNG image as uint16
img = imread('../../data/scream/scream_rgb.png');  % Replace with your image path

% Ensure the image is uint16 (If it's not, convert it)
if ~isa(img, 'uint16')
    img = uint16(img);  % Convert to uint16 if not already
end

% Display the image
figure;
imshow(img);
title('Select an area to modify');

% Let the user select a rectangular region using drawrectangle
h = drawrectangle;
wait(h);  % Wait for the user to select the rectangle

% Get the position of the selected rectangle
position = h.Position;  % position = [x, y, width, height]

% Extract the region of interest (ROI) based on the selected rectangle
x1 = round(position(1));
y1 = round(position(2));
x2 = round(position(1) + position(3));
y2 = round(position(2) + position(4));

% Extract the selected region from the image
selected_region = img(y1:y2, x1:x2, :);

% Convert the selected region to HSV color space to modify saturation and value
selected_region_hsv = rgb2hsv(double(selected_region) / 65535);  % Normalize to [0, 1] before conversion

% Reduce saturation by a factor (e.g., 0.5 for 50% reduction)
selected_region_hsv(:, :, 2) = selected_region_hsv(:, :, 2) * 0.5;  % Saturation

% Reduce contrast by modifying the value channel (0.7 factor for 30% reduction)
selected_region_hsv(:, :, 3) = selected_region_hsv(:, :, 3) * 0.7;  % Value/brightness

% Ensure values are clamped between [0, 1]
selected_region_hsv(:, :, 2) = max(0, min(1, selected_region_hsv(:, :, 2)));  % Saturation (S)
selected_region_hsv(:, :, 3) = max(0, min(1, selected_region_hsv(:, :, 3)));  % Value (V)

% Convert back to RGB
modified_region_rgb = hsv2rgb(selected_region_hsv);

% Convert the modified region back to uint16 (scaling back to [0, 65535])
modified_region_rgb = uint16(modified_region_rgb * 65535);

% Create a linear gradient mask for smooth blending
[rows, cols, ~] = size(modified_region_rgb);

% Create a mask that gradually transitions across the width or height of the region
gradient_mask = linspace(0, 1, cols);  % Linear gradient along the width
gradient_mask = repmat(gradient_mask, rows, 1);  % Extend to the full height

% Now, blend the original and modified regions using the gradient mask
blended_region = uint16(zeros(size(modified_region_rgb)));

for c = 1:3  % For each color channel (R, G, B)
    % Applying the linear mask along the width
    blended_region(:, :, c) = uint16(gradient_mask .* double(modified_region_rgb(:, :, c)) + ...
                                     (1 - gradient_mask) .* double(selected_region(:, :, c)));
end

% Ensure blended region stays within valid uint16 range
blended_region = min(max(blended_region, 0), 65535);  % Clamping to [0, 65535]

% Replace the modified and blended region in the original image
img(y1:y2, x1:x2, :) = blended_region;

% Display the modified image
figure;
imshow(img, []);  % Display with appropriate scaling
title('Image with Reduced Saturation and Contrast in Selected Region');

% Optionally, save the modified image
imwrite(img, 'modified_image_blend_stripe.png');  % Save the modified image
