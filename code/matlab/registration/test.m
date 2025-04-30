% Create a simple gradient image
I = repmat(linspace(0, 1, 256), 256, 1);

% Apply imadjust with default values (no input range specified)
I_adjusted = imadjust(I);

% Display the original and adjusted images
subplot(1, 2, 1);
imshow(I, []);
title('Original Image');

subplot(1, 2, 2);
imshow(I_adjusted, []);
title('Image after imadjust');
