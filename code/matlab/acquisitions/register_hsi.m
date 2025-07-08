%% Clear and load your two already-stitched hypercubes (moving & fixed)
clear; close all;

% ——— Filepaths —————————————————————————————————————————————
hdrMov = "/home/oem/eliza/data/processed/reflectance/after/yoda_reflectance_after_full.hdr";
hdrFix = "/home/oem/eliza/data/processed/reflectance/before/yoda_reflectance_full_before.hdr";

% ——— Load hypercubes —————————————————————————————————————
hcMov = hypercube(hdrMov);
hcFix = hypercube(hdrFix);

% ——— Band mask (380–780 nm) ———————————————————————————————
wlMask = hcMov.Wavelength >= 380 & hcMov.Wavelength <= 780;
bands  = find(wlMask);
wl     = hcMov.Wavelength(bands);
nBands = numel(bands);

% ——— Extract & gather the data cubes ————————————————————————
movCube = gather(hcMov.DataCube(:,:,bands));   % [rows × cols × nBands]
fixCube = gather(hcFix.DataCube(:,:,bands));

%% 1) Draw ROI on band 100 of the MOVING cube
refBandIndex = 100;    % use band 100
threshold    = 1;      % clip level

% Prepare the band for display
bandImg = double(movCube(:,:,refBandIndex));
bandImg(bandImg > threshold) = threshold;
bandImg = (bandImg - min(bandImg(:))) / (max(bandImg(:)) - min(bandImg(:)));
bandImgAdj = imadjust(bandImg);

% Display and draw ROI with drawrectangle
figure; imshow(bandImgAdj, []);
axis image off;
title('Draw ROI on MOVING cube (band 100)','FontSize',14);
hRec = drawrectangle('StripeColor','r');
wait(hRec);
pos = round(hRec.Position);  % [xmin, ymin, width, height]

% Crop MOVING cube across all bands
movCrop = movCube( ...
    pos(2):(pos(2)+pos(4)-1), ...
    pos(1):(pos(1)+pos(3)-1), ...
    : );

% FIXED cube remains full-frame
fixCrop = fixCube;

%% 2) Prepare the reference band for registration
movRef = double(movCrop(:,:,refBandIndex));
fixRef = double(fixCrop(:,:,refBandIndex));

% Clip, normalize, contrast-stretch
movRef(movRef > threshold) = threshold;
fixRef(fixRef > threshold) = threshold;
movRef = (movRef - min(movRef(:))) / (max(movRef(:)) - min(movRef(:)));
fixRef = (fixRef - min(fixRef(:))) / (max(fixRef(:)) - min(fixRef(:)));
movRefAdj = imadjust(movRef);
fixRefAdj = imadjust(fixRef);

%% 3) Detect & match features, estimate affine transform
ptsMov  = detectSURFFeatures(movRefAdj);
ptsFix  = detectSURFFeatures(fixRefAdj);
[featMov, validMov] = extractFeatures(movRefAdj, ptsMov);
[featFix, validFix] = extractFeatures(fixRefAdj, ptsFix);
idxPairs = matchFeatures(featMov, featFix, 'MaxRatio',0.8, 'Unique',true);

matchedMov = validMov(idxPairs(:,1));
matchedFix = validFix(idxPairs(:,2));
[tformObj, inlierIdx] = estimateGeometricTransform2D( ...
    matchedMov, matchedFix, 'affine', ...
    'MaxDistance',1.8, 'Confidence',99, 'MaxNumTrials',2000);
%%
% Visualize inlier matches
figure;
showMatchedFeatures(fixRefAdj, movRefAdj, matchedFix(inlierIdx), matchedMov(inlierIdx), 'montage');
title('Inlier Matches','FontSize',14);

%% 4) Define outputView = the fixed-image frame (registration only)
Rfixed     = imref2d(size(fixRef));   % match fixed band size
outputView = Rfixed;

%% 5) Warp every band of the MOVING cube into the fixed frame
registeredMov = zeros(size(fixCube), 'like', fixCube);  % [rows × cols × nBands]

for b = 1:nBands
    wMov = imwarp(movCrop(:,:,b), tformObj, 'OutputView', outputView);
    registeredMov(:,:,b) = wMov;
end

%% 6) False‐RGB of FIXED vs. REGISTERED MOVING
fixCube1 = fixCube;
fixCube1(fixCube1>1) = 1;
registeredMov1 = registeredMov;
registeredMov1(registeredMov1>1)=1;
hcFixed   = hypercube(fixCube1,      wl);
hcRegMov  = hypercube(registeredMov1, wl);

rgbFixed  = colorize(hcFixed,  'Method','RGB','ContrastStretching',true);
rgbRegMov = colorize(hcRegMov, 'Method','RGB','ContrastStretching',true);
%%
figure;
imshowpair(rgbFixed, rgbRegMov, 'montage');
title('Fixed (left) vs. Registered Moving (right)','FontSize',16);

%% 7) Overlay to evaluate registration quality
% a) Blend overlay
figure;
imshowpair(rgbFixed, rgbRegMov);
title('Overlay false RGB','FontSize',16);
%%
% b) False‐color overlay (fixed=red, registered=green)
overlayFC = imfuse(rgbFixed, rgbRegMov, ...
    'falsecolor', 'Scaling','independent', 'ColorChannels',[1 2 0]);
figure;
imshow(overlayFC);
title('Overlay (falsecolor): Red=Fixed, Green=Registered','FontSize',16);
%
%% ——— 8) Overlay fixed vs. warped-moving on the reference band ——————
% warp the full moving-band 100 into the fixed frame
wMovRef = imwarp(movRefAdj, tformObj, 'OutputView', outputView );

% grab the fixed band
fixRefFull = fixRefAdj;

% a) Blend overlay
figure;
imshowpair(fixRefFull, wMovRef);
colormap gray;
title('Overlay Fixed vs. Warped Moving — Band 100','FontSize',14);
%%
% b) False-color overlay (red = fixed, green = warped-moving)
fc = imfuse(fixRefFull, wMovRef, ...
    'falsecolor','Scaling','independent','ColorChannels',[1 2 0]);
figure;
imshow(fc);
title('Overlay (falsecolor) of Fixed (red) vs. Warped Moving (green) — Band 100','FontSize',14);

%% save the registered cube


outFolder = "/home/oem/eliza/data/processed/reflectance/registered";
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end
imgPath = fullfile(outFolder, "yoda_reflectance_after_reg.img");
hdrPath = fullfile(outFolder, "yoda_reflectance_after_reg.hdr");

% b) Clip to [0,1] if needed & convert to single precision
regCube = single(registeredMov);

[rows, cols, bands] = size(regCube);

% c) Write the binary cube (.img) in BSQ order
multibandwrite( ...
    regCube, ...           % data array
    imgPath, ...           % output .img path
    'bsq', ...             % band-sequential
    'Precision','single', ...  % 32-bit float
    'Offset',0 ...            % no header offset
);

% d) Write the ENVI header (.hdr) with wavelength metadata
fid = fopen(hdrPath, 'w');
fprintf(fid, 'ENVI\n');
fprintf(fid, 'samples = %d\n', cols);
fprintf(fid, 'lines   = %d\n', rows);
fprintf(fid, 'bands   = %d\n', bands);
fprintf(fid, 'header offset = 0\n');
fprintf(fid, 'file type = ENVI Standard\n');
fprintf(fid, 'data type = 4\n');          % 4 == 32-bit float
fprintf(fid, 'interleave = bsq\n');
fprintf(fid, 'byte order = 0\n');         % little-endian
fprintf(fid, 'wavelength units = nm\n');

% write the wavelength list
fprintf(fid, 'wavelength = {');
for k = 1:bands
    fprintf(fid, '%g', wl(k));
    if k < bands
        fprintf(fid, ',');
    end
end
fprintf(fid, '}\n');

fclose(fid);

% e) (Optional) reload as a hypercube object to confirm
hcReg = hypercube(hdrPath);

% f) Quick false-RGB check
rgbReg = colorize(hcReg, 'Method','RGB', 'ContrastStretching', true);
figure; imshow(rgbReg);
title('Registered Moving Cube — False-RGB','FontSize',16);


%% Extract specularity mask

% 1) Load your registered hypercube and pull the data
clear;
hc   = hypercube("/home/oem/eliza/data/reflectance/registered/cubes/cactus_reflectance_full_before.hdr");
data_cube  = gather(hc.DataCube);   % [rows × cols × bands]
wl = hc.Wavelength;

%% 
minRefl = min(data_cube, [], 3);
%
% 1) Plot using histogram with explicit bin edges
edges = 0:0.01:1;
h = histogram(minRefl(:), 'BinEdges', edges, ...
              'FaceColor','b', 'EdgeColor','none');

% Make the axis labels big and bold
xlabel('Min reflectance across bands',  ...
       'FontSize',14, 'FontWeight','bold');
ylabel('Pixel count',                      ...
       'FontSize',14, 'FontWeight','bold');

% Grab the current axes and bump up the tick‐label size & weight
ax = gca;
ax.FontSize   = 20;        % larger tick labels
ax.FontWeight = 'bold';    % bold tick labels
%%
%
for T = [0.7 0.8 0.9]
    m = minRefl >= T;
    fprintf('Threshold = %.2f: %.2f%% of pixels flagged\n', ...
        T, 100*nnz(m)/numel(m));
end


%%
% 2) Choose your specularity threshold
thresh = 0.55;

% 3) Build the 2-D mask: true where every band ≥ thresh
spec_mask = all(data_cube >= thresh, 3);   % [rows × cols] logical

% 4) (Optional) Clean it up by removing tiny speckles
%    specMask2D = bwareaopen(specMask2D, 50);

% 5) Visualize it
figure;
imshow(spec_mask);
title(sprintf('Specularity Mask Cactus Before Ageing', thresh), 'FontSize',14);
colormap gray; axis image off;
%%
data_cube(data_cube>1) = 1;
outMat = "/home/oem/eliza/data/reflectance/registered/cactus_reflectance_before.mat";
save(outMat, 'data_cube', 'wl', 'spec_mask', '-v7.3');


%% Fast sRGB computation from reflectance hypercube
fprintf('--- Computing sRGB from hyperspectral data ---\n');

% Load CIE 1931 2-deg CMFs and D50 illuminant
cmf_path = '../../../data/CIE2degCMFs_full.csv';  % Wavelength in col 1, x y z in 2:4
D50_path = '../../../data/CIE_D50.txt';           % Wavelength in col 1, SPD in col 2

cmf_data = importdata(cmf_path);
D50_data = importdata(D50_path);

cmf_wl = cmf_data(:,1);
cmf_vals = cmf_data(:,2:4);  % x̄, ȳ, z̄
D50_wl = D50_data(:,1);
D50_vals = D50_data(:,2);    % Spectral power distribution

% Interpolate CMFs and illuminant to match your wavelengths
cmf_interp = interp1(cmf_wl, cmf_vals, wl, 'linear', 'extrap');     % [#bands × 3]
D50_interp = interp1(D50_wl, D50_vals, wl, 'linear', 'extrap');     % [#bands × 1]

% Reshape hyperspectral cube to 2D [N_pixels × bands]
[H, W, B] = size(data_cube);
refl_2d = reshape(data_cube, [], B);

% Compute XYZ tristimulus values
XYZ = ref2xyz(D50_interp(:), cmf_interp, refl_2d);  % [N_pixels × 3]

% Convert XYZ to sRGB (D65 adaptation assumed)
RGB = xyz2rgb(XYZ./100, 'ColorSpace', 'srgb', ...
                    'WhitePoint', 'd50', ...
                    'OutputType', 'double');
% Reshape back to image
RGB_img = reshape(RGB, H, W, 3);
RGB_img = min(max(RGB_img, 0), 1);  % Clip for display

% Show the sRGB image
figure;
imshow(RGB_img);
title('sRGB Image from Hyperspectral Cube');

