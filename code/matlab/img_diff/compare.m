
clc; clear; close all;
roof = double(intmax('uint16'));

rng(10);

cubeFile = "/Volumes/School/Thesis/Data/Scream/1971/Ekta1971_restored.hdr";
refFile = "/Volumes/School/Thesis/Data/Scream/present/inpainted_scream.hdr";

ill = importdata('../../../data/CIE_D65.txt'); % Illuminant
CMFs = importdata('../../../data/CIE2degCMFs_1931.txt'); % Color matching functions

hcube = hypercube(cubeFile);
inCUBE = hcube.DataCube;
bands = hcube.Wavelength;
[m, n, bd] = size(inCUBE);

hcube_ref = hypercube(refFile);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;
[mr, nr, bd_ref_count] = size(refCUBE);

% cropping the reference cube to 380 - 780 nm
valid_idx = find(bands_ref >= 380 & bands_ref <= 780);

refCUBE = refCUBE(:, :, valid_idx);
bands_ref = bands_ref(valid_idx);
bd_ref_count = length(valid_idx);


% Flatten the spectral data
lincube = double(reshape(inCUBE, [], bd));
lincube_ref = reshape(refCUBE, [], bd_ref_count);


%% Interpolate illuminant and CMFs for the input cube
illIP = interp1(ill(:,1), ill(:,2), bands, 'spline');
CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,3), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,4), bands, 'spline')];
sp_tristREF = CMFsIP .* illIP;
xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);

% For reference cube:
illIP_ref = interp1(ill(:,1), ill(:,2), bands_ref, 'spline');
CMFsIP_ref = [interp1(CMFs(:,1), CMFs(:,2), bands_ref, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,3), bands_ref, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,4), bands_ref, 'spline')];
sp_tristREF_ref = CMFsIP_ref .* illIP_ref;
xyz_ref = (lincube_ref * sp_tristREF_ref) ./ sum(sp_tristREF_ref(:,2), 1);


%% ---------------- Convert XYZ to RGB ----------------
% Convert using the same color space as your display.
rgb_input = xyz2rgb(xyz_input, 'ColorSpace', 'prophoto-rgb');
rgb_ref   = xyz2rgb(xyz_ref,   'ColorSpace', 'prophoto-rgb');

% Reshape into 10x14 for patch indexing
rgb_input_2d = reshape(rgb_input, m, n, 3);
rgb_ref_2d   = reshape(rgb_ref,   mr, nr, 3);
rgb_input_uint16 = im2uint16(mat2gray(rgb_input_2d));
rgb_ref_uint16   = im2uint16(mat2gray(rgb_ref_2d));

%%
% Display
% imshow(rgb_input_uint16);
% figure;
% imshow(rgb_ref_uint16);

% imwrite(rgb_input_uint16, 'results/rgb_input.png');
% imwrite(rgb_ref_uint16, 'results/rgb_ref.png');

%% Cropping and Registration

% cropped_input = imread("results/rgb_ref.png");
% cropped_ref = imread("results/rgb_ref.png");


%%
figure('Name','Select ROI for Input Image','NumberTitle','off');
imshow(rgb_input_uint16);
title('Draw a rectangle to crop the input image');
hInput = drawrectangle();          % interactive ROI
rectInput = round(hInput.Position);     % [x, y, width, height]
input_cropped = imcrop(rgb_input_uint16, rectInput);
close;

% 2) Manually Crop the Reference Image
figure('Name','Select ROI for Reference Image','NumberTitle','off');
imshow(rgb_ref_uint16);
title('Draw a rectangle to crop the reference image');
hRef = drawrectangle();          % interactive ROI
rectRef = round(hRef.Position);     % [x, y, width, height]
ref_cropped = imcrop(rgb_ref_uint16, rectRef);
close;

%%
[h1, w1, ~] = size(input_cropped);
[h2, w2, ~] = size(ref_cropped);

target_size = [min(h1, h2), min(w1, w2)];
input_resized = imresize(input_cropped, target_size);
ref_resized   = imresize(ref_cropped,   target_size);

[registered, tform] = register_images(input_resized, ref_resized);
% --- Remove black padding caused by imwarp ---
mask = sum(registered, 3) > 0; % Find non-black areas
props = regionprops(mask, 'BoundingBox');
if ~isempty(props)
    bbox = round(props(1).BoundingBox);
    registered = imcrop(registered, bbox);         % Crop registered
    ref_resized = imcrop(ref_resized, bbox);       % Optionally crop reference too
end

imwrite(registered, 'results/rgb_input_registered.png');

%% Using Registration Estimator

moving_img = rgb2gray(input_cropped);
fixed_img = rgb2gray(ref_cropped);
registrationEstimator(moving_img, fixed_img);
f = warndlg(['Load images from workspace. Export registered image ' ...
    'with the default name. When finished press OK.'], 'Registration');
waitfor(f);

%%

inputIMGreg = uint16(zeros(size(fixed_img,1),size(fixed_img,2),3));
fixedRef = imref2d(movingReg.SpatialRefObj.ImageSize);
for i=1:3
    inputIMGreg(:,:,i) = imwarp(input_cropped(:,:,i),...
        movingReg.Transformation,'OutputView',fixedRef);
end

%%
diff_image = imabsdiff(ref_cropped, inputIMGreg);
imshow(diff_image, []);


%% The black pixel at the edge of the registered image are cropped out
g = rgb2gray(inputIMGreg);
g = [g(1,:) g(end,:) g(:,1)' g(:,end)'];
while any(g==0)
    inputIMGreg = inputIMGreg(2:(end-1),2:(end-1),:);
    g = rgb2gray(inputIMGreg);
    g = [g(1,:) g(end,:) g(:,1)' g(:,end)'];
end
crop = (size(ref_cropped,1)-size(inputIMGreg,1))/2;
refIMGcrop = ref_cropped(crop+1:end-crop,crop+1:end-crop,:);

diff_image = imabsdiff(refIMGcrop, inputIMGreg);
imshow(diff_image, []);

%% Registration

% First registration to estimate scaling factor
% [registered1, tform1, scaling_factor] = register_images(cropped_input, cropped_ref);
% fprintf('Estimated pixel spacing ratio (reference/input): %.4f\n', scaling_factor);
% 
% % Resize the input image based on the scaling factor
% resized_input = imresize(cropped_input, scaling_factor);
% 
% % Second registration after resizing
% [registered2, tform2, ~] = register_images(resized_input, cropped_ref);

%%




function [registered, tform] = register_images(input_image, ref_image)
    % Convert to grayscale
    gray1 = im2single(rgb2gray(input_image));
    gray2 = im2single(rgb2gray(ref_image));

    % 1) Detect & extract SURF features
    pts1 = detectSURFFeatures(gray1, 'MetricThreshold', 500);      % Adjust threshold as needed :contentReference[oaicite:0]{index=0}
    pts2 = detectSURFFeatures(gray2, 'MetricThreshold', 500);      % :contentReference[oaicite:1]{index=1}
    [f1, vpts1] = extractFeatures(gray1, pts1);
    [f2, vpts2] = extractFeatures(gray2, pts2);

    % 2) Match features
    indexPairs = matchFeatures(f1, f2, ...
        'MatchThreshold', 50, ...      % Tighter threshold for SURF
        'MaxRatio', 0.7, ...
        'Unique', true);
    matched1 = vpts1(indexPairs(:,1));
    matched2 = vpts2(indexPairs(:,2));

    % 3) Estimate projective transform from inliers
    [tform, inlierIdx] = estimateGeometricTransform2D( ...
        matched1.Location, matched2.Location, 'projective', ...
        'MaxDistance', 4, ...
        'Confidence', 99.9, ...
        'MaxNumTrials', 3000);

    % 4) Warp the cropped input to the reference frame
    outputView = imref2d(size(ref_image(:,:,1)));
    registered = imwarp(input_image, tform, 'OutputView', outputView);

    % 5) Compute scaling factor (for optional further resizing)
    loc1 = matched1.Location;
    loc2 = matched2.Location;
    d_input = sqrt(sum((loc1 - mean(loc1)).^2, 2));
    d_ref   = sqrt(sum((loc2 - mean(loc2)).^2, 2));

    % (Optional) Display results
    inlierPts1 = matched1(inlierIdx);
    inlierPts2 = matched2(inlierIdx);
    figure('Name','SURF Registration','NumberTitle','off','Position',[100 100 1200 800]);
    subplot(2,3,1), imshow(input_image), title('1: Input');
    subplot(2,3,2), imshow(ref_image),   title('2: Reference');
    subplot(2,3,3), showMatchedFeatures(gray1, gray2, inlierPts1, inlierPts2, 'montage');
                     title('3: Inlier SURF Matches');
    subplot(2,3,4), imshowpair(ref_image, registered);
                     title('4: Registered vs. Reference');
    diff_image = imabsdiff(ref_image, registered);
    subplot(2,3,5), imshow(diff_image, []), title('5: Abs Difference');
end
