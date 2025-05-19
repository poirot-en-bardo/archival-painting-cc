ms_img_folder = "../data/ChinaGirl/7_set";
led_path = '../data/ChinaGirl/CVCL7bands.txt';
cmf_path = '../../../data/CIE2degCMFs_full.csv';

% Get all TIFF filenames in the folder:
files = dir(fullfile(ms_img_folder, '*.tif'));
assert(numel(files)==7, 'Expected 7 TIFFs in %s, found %d.', ms_img_folder, numel(files));

% Read the first image to grab size + class:
firstImg = imread(fullfile(ms_img_folder, files(1).name));
[H, W, C] = size(firstImg);
assert(C==1, 'Expected single‐band (grayscale) TIFFs.');

% Preallocate cube (H×W×7) of the same class:
cube = zeros(H, W, numel(files), class(firstImg));

% Loop through and fill:
for k = 1:numel(files)
    img = imread(fullfile(ms_img_folder, files(k).name));
    assert(isequal(size(img), [H W]), ...
           'Image %s is a different size.', files(k).name);
    cube(:,:,k) = img;
end

%% checking the cube
% imshow(cube(:,:,3), []);

LEDset = readmatrix(led_path, 'Delimiter','\t');
wl_led = LEDset(:,1);       

fullCMFs = importdata(cmf_path);

wl_cmf = fullCMFs(:, 1);
idx = wl_cmf >= 380 & wl_cmf <= 780;
wl_cmf = wl_cmf(idx);      
cmfs_vis     = fullCMFs(idx, :);  
CMFs = (cmfs_vis(:,2:end)'*LEDset(:,2:end)./sum(LEDset(:,2:end)))';


[r, c, ~] = size(cube);
lincube_capture = reshape(double(cube)/roof,[],bandN);
xyz_MS_full= (lincube_capture * CMFs)./sum(CMFs(:,2),1);

rgbImg_MS_full = xyz2rgb(xyz_MS_full,'ColorSpace','prophoto-rgb');
rgbImg_MS_full  = reshape(rgbImg_MS_full, r, c, 3);  

% white-balance by dividing by the white patch
% whiteIdx = 7;        % row 2, col 1 → patch #7
% rgbMSwb_actual = rgbMS_actual ./ rgbMS_actual(whiteIdx, :);
p = 7;
x = rects(p,1);  
y = rects(p,2);  
w = rects(p,3);  
h = rects(p,4);  

% turn into integer pixel indices, clamp:
xs = max(1,floor(x)) : min(size(rgbImg_MS_full,2), ceil(x+w));
ys = max(1,floor(y)) : min(size(rgbImg_MS_full,1), ceil(y+h));

% extract that patch across all 3 channels and average:
patch_rgb     = rgbImg_MS_full(ys, xs, :);            % H×W×3
white_patch   = squeeze( mean(mean(patch_rgb,1),2) );  % 3×1
%%
% now white‐balance the whole image:
rgbWB = rgbImg_MS_full ./ reshape(white_patch, 1,1,3);
rgbWB = max(min(rgbWB,1),0);   % clamp
rgb16   = im2uint16(rgbWB);  

outFile = fullfile(outDir,'capture_act_full.tif');
% saveProPhotoTIFF(rgb16, outFile);

% display
figure; imshow(rgb16);
title('Actual MS Capture (RGB) - full');