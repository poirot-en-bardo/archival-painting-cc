%% Clear and load your two already-stitched hypercubes (moving & fixed)
clear; close all;

% ——— Filepaths —————————————————————————————————————————————
matPath   = "/home/oem/eliza/data/film_scans/white_balanced/yoda_halogen_fuji_underexp_balanced.mat";
hdrFix    = "/home/oem/eliza/data/processed/reflectance/registered/yoda_reflectance_full_before.hdr";
refBandHS = 100;   % Band from HSI (fixed)
refBandMS = 9;     % Band from MS (moving)

% ——— Load rcubes —————————————————————————————————————
S       = load(matPath, 'wb_cube', 'wl');
movCube = flip(S.wb_cube, 1);             % [rows × cols × M]
wlMS    = S.wl;

hcFix   = hypercube(hdrFix);
fixCube = gather(hcFix.DataCube);         % [rows × cols × H]
wlFix   = hcFix.Wavelength;

%% 1) Draw ROI on band 100 of the MOVING cube
threshold = 1;

% Prepare the band for display
bandImg = double(movCube(:,:,refBandMS));
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
movRef = double(movCrop(:,:,refBandMS));
fixRef = double(fixCrop(:,:,refBandHS));

% Clip, normalize, contrast-stretch
movRef(movRef > threshold) = threshold;
fixRef(fixRef > threshold) = threshold;
movRef = (movRef - min(movRef(:))) / (max(movRef(:)) - min(movRef(:)));
fixRef = (fixRef - min(fixRef(:))) / (max(fixRef(:)) - min(fixRef(:)));
%%
movRefAdj = imadjust(movRef, [], [], 0.6);
fixRefAdj = imadjust(fixRef);

%% 3) Detect & match features using SIFT
%% Detect and combine features from multiple detectors (fixed)
% --- MOVING IMAGE ---
% --- MOVING IMAGE ---
ptsMov = detectSIFTFeatures(movRefAdj, 'ContrastThreshold', 0.001);  % lower threshold = more features

% --- FIXED IMAGE ---
ptsFix = detectSIFTFeatures(fixRefAdj, 'ContrastThreshold', 0.001);


% --- Feature Extraction ---
[featMov, validMov] = extractFeatures(movRefAdj, ptsMov);
[featFix, validFix] = extractFeatures(fixRefAdj, ptsFix);

%
% --- Feature Matching ---
idxPairs = matchFeatures(featMov, featFix, ...
    'MaxRatio', 0.8, ...
    'MatchThreshold', 100, ...
    'Unique', true);

matchedMov = validMov(idxPairs(:,1));
matchedFix = validFix(idxPairs(:,2));

% --- Estimate Transform ---
[tformObj, inlierIdx] = estimateGeometricTransform2D( ...
    matchedMov, matchedFix, 'projective', ...
    'MaxDistance', 20, 'Confidence', 97, 'MaxNumTrials', 10000);
%%
% --- Visualize Matches ---
figure;
showMatchedFeatures(fixRefAdj, movRefAdj, ...
    matchedFix(inlierIdx), matchedMov(inlierIdx), 'montage');
title('Inlier Matches (Combined Detectors)', 'FontSize', 14);



%% 4) Define outputView = the fixed-image frame (registration only)
Rfixed     = imref2d(size(fixRef));
outputView = Rfixed;

%% 5) Warp every band of the MOVING cube into the fixed frame
registeredMov = zeros(size(fixCube), 'like', fixCube);
[~, ~, nBandsMov]  = size(movCrop);
[rowsFix, colsFix] = size(fixRef);
registeredCube = zeros(rowsFix, colsFix, nBandsMov, 'like', movCrop);

for b = 1:nBandsMov
    wMov = imwarp(movCrop(:,:,b), tformObj, 'OutputView', outputView);
    registeredCube(:,:,b) = wMov;
end

%% 6) Overlay fixed vs. warped-moving on the reference band
wMovRef = imwarp(movRefAdj, tformObj, 'OutputView', outputView);
fixRefFull = fixRefAdj;

% a) Blend overlay
figure;
imshowpair(fixRefFull, wMovRef);
colormap gray;
title('Overlay Fixed vs. Warped Moving — Band 100 (SIFT)','FontSize',14);

