rootFolder = '/home/oem/eliza/data/film_scans/underexp/';
outputDir  = '/home/oem/eliza/mac-shared/film_scans';
led_path = '../../../data/film/CVCL10bands.txt';
cmf_path = '../../../data/CIE2degCMFs_full.csv';
D50_path = '../../../data/CIE_D50.txt';

D50_SPD = importdata(D50_path);  
fullCMFs = importdata(cmf_path);
LEDset = readmatrix(led_path, 'Delimiter','\t');
wl_led = LEDset(:,1);       

roof = double(intmax('uint16'));
bandN = size(LEDset, 2) - 1;  

% List all immediate subfolders of rootFolder (excluding . and ..)
allChildren = dir(rootFolder);
allChildren = allChildren([allChildren.isdir]);
allChildren = allChildren(~ismember({allChildren.name}, {'.', '..'}));

% for all subfolders
for i = 1:numel(allChildren)
    childName = allChildren(i).name;
    flatfieldedPath = fullfile(rootFolder, childName, 'FlatFielded');
    
    if ~exist(flatfieldedPath, 'dir')
        warning('No FlatFielded folder in %s, skipping.', childName);
        continue;
    end
    
    ms_img_folder = flatfieldedPath; 

    %% Load the multispectral cube

    % Get all TIFF filenames in the folder:
    files = dir(fullfile(ms_img_folder, '*.tif'));

    % Read the first image to grab size + class:
    firstImg = imread(fullfile(ms_img_folder, files(1).name));
    [H, W, C] = size(firstImg);

    % Preallocate cube (H×W×7) of the same class:
    cube = zeros(H, W, numel(files), class(firstImg));

    % Loop through and fill:
    for k = 1:numel(files)
        img = imread(fullfile(ms_img_folder, files(k).name));
        cube(:,:,k) = img;
    end

    %% create the RGB rendering

    d50_spd = interp1(D50_SPD(:,1), D50_SPD(:,2), wl_led, 'linear', 'extrap');
    cmf_xyz = interp1(fullCMFs(:,1), fullCMFs(:,2:4), wl_led, 'linear', 'extrap');

    d50_band = zeros(bandN, 1);    
    cmf_band = zeros(bandN, 3);     

    % weight
    for k = 1:bandN
        led_spd = LEDset(:, k+1);    
        norm_led = led_spd / sum(led_spd);  

        d50_band(k)      = sum(d50_spd   .* norm_led);
        cmf_band(k, 1)   = sum(cmf_xyz(:,1) .* norm_led); 
        cmf_band(k, 2)   = sum(cmf_xyz(:,2) .* norm_led); 
        cmf_band(k, 3)   = sum(cmf_xyz(:,3) .* norm_led); 
    end

    cube_dbl = double(cube) / double(intmax('uint16'));  
    cube_illum = cube_dbl .* reshape(d50_band, 1, 1, bandN);

    [nRows, nCols, ~] = size(cube_illum);
    cube_flat = reshape(cube_illum, [], bandN);

    XYZ = cube_flat * cmf_band;         

    k = sum(d50_band .* cmf_band(:,2));
    XYZ = XYZ / k;
    XYZimg   = reshape(XYZ, H, W, 3);

    RGB_pPhoto = xyz2rgb(XYZimg, 'ColorSpace', 'prophoto-rgb', 'WhitePoint','d50');
    RGB_linear = xyz2rgb(XYZimg, 'ColorSpace', 'linear-rgb', 'WhitePoint','d65');


    %% White balancing based on CC white patch
    figure; imshow(RGB_pPhoto);
    title('Draw rectangle around white patch and double-click');
    hW = drawrectangle();            
    wait(hW);                 % block until double-click
    whitePos = hW.Position;

    xw    = max(1, floor(whitePos(1)));
    yw    = max(1, floor(whitePos(2)));
    ww    = floor(whitePos(3));
    hw    = floor(whitePos(4));
    xw_end = min(size(RGB_pPhoto,2), xw + ww  - 1);
    yw_end = min(size(RGB_pPhoto,1), yw + hw  - 1);

    patch_rgb   = RGB_pPhoto(yw:yw_end, xw:xw_end, :); 
    white_patch = squeeze( mean(mean(patch_rgb,1),2) );    % average

    % white-balance the whole image
    rgbWB = RGB_pPhoto ./ reshape(white_patch,1,1,3);
    rgbWB = max(min(rgbWB,1),0);   % [0,1]
    rgb16 = im2uint16(rgbWB);

    %% cropping the image
    figure; imshow(rgb16);
    title('Draw rectangle to crop and double-click');
    hC = drawrectangle();            
    wait(hC);                 
    cropPos = hC.Position;    

    % clamp & round
    xc    = max(1, floor(cropPos(1)));
    yc    = max(1, floor(cropPos(2)));
    wc    = floor(cropPos(3));
    hc2   = floor(cropPos(4));
    xc_end = min(size(rgb16,2), xc + wc - 1);
    yc_end = min(size(rgb16,1), yc + hc2 - 1);

    rgbCropped = rgb16(yc:yc_end, xc:xc_end, :);
    rgbCropped_linear = im2uint16(RGB_linear(yc:yc_end, xc:xc_end, :));


    figure; imshow(rgbCropped);
    title('Cropped RGB Image');

    %% Save the cropped image

    % Construct output file name based on child folder
    outFile = fullfile(outputDir, [childName '_d50.png']);

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    imwrite(rgbCropped, outFile);
    imwrite(rgbCropped_linear, fullfile(outputDir, [childName 'linear_rgb.png']));
    fprintf('Saved: %s\n', outFile);
end
