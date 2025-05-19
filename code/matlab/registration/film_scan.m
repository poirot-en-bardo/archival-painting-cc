ms_img_folder = "../data/ChinaGirl/7_set";
sensor_path = '../data/ChinaGirl/IMX455_QE.txt'; 
tsr_path = '../data/ChinaGirl/set2_black.csv';
led_path = '../data/ChinaGirl/CVCL7bands.txt';
cmf_path = '../data/CIE2degCMFs_full.csv';

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