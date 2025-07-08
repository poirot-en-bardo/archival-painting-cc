clear; close all;

% --- Filepaths ---
matPath   = "/home/oem/eliza/data/film_scans/white_balanced/cactus_halogen_kodak_exp0_balanced.mat";
hdrFix    = "/home/oem/eliza/data/processed/reflectance/registered/cactus_reflectance_full_before.hdr";
refBandHS = 100;   % Band from HSI (fixed)
refBandMS = 9;     % Band from MSI (moving)

% --- Load cubes ---
S       = load(matPath, 'wb_cube', 'wl');
movCube = flip(S.wb_cube, 1);             
wlMS    = S.wl;
hcFix   = hypercube(hdrFix);
fixCube = gather(hcFix.DataCube);         
wlFix   = hcFix.Wavelength;

% --- ROI on MOVING cube ---
threshold = 1;
bandImg = double(movCube(:,:,refBandMS));
bandImg(bandImg > threshold) = threshold;
bandImg = (bandImg - min(bandImg(:))) / (max(bandImg(:)) - min(bandImg(:)));
bandImgAdj = imadjust(bandImg);

figure; imshow(bandImgAdj, []);
axis image off;
title('Draw ROI on MOVING cube (band 9)','FontSize',14);
hRec = drawrectangle('StripeColor','r');
wait(hRec);
pos = round(hRec.Position);

movCrop = movCube( ...
    pos(2):(pos(2)+pos(4)-1), ...
    pos(1):(pos(1)+pos(3)-1), ...
    : );
fixCrop = fixCube;

% --- Prepare registration bands ---
movRef = double(movCrop(:,:,refBandMS));
fixRef = double(fixCrop(:,:,refBandHS));
movRef(movRef > threshold) = threshold;
fixRef(fixRef > threshold) = threshold;
movRef = (movRef - min(movRef(:))) / (max(movRef(:)) - min(movRef(:)));
fixRef = (fixRef - min(fixRef(:))) / (max(fixRef(:)) - min(fixRef(:)));
movRefAdj = imadjust(movRef, [], [], 0.8);
fixRefAdj = imadjust(fixRef);

% --- Save images and launch Registration Estimator ---
imwrite(im2uint8(fixRefAdj), 'fixRefAdj.png');
imwrite(im2uint8(movRefAdj), 'movRefAdj.png');
disp("Launching Registration Estimator. Load fixRefAdj.png as Fixed, movRefAdj.png as Moving.");
registrationEstimator;

% Wait for user to finish registration and export 'tform'
disp("Waiting for you to finish registration in the GUI and export the 'tform' variable.");
disp("Press any key to continue once 'tform' is available in the workspace.");
pause;

% --- Check for 'tform' ---
if ~exist('tform', 'var')
    error("'tform' variable not found. Make sure you exported it from Registration Estimator.");
end

% --- Use exported Rfixed or default to fixed image size ---
if exist('Rfixed', 'var')
    outputView = Rfixed;
else
    outputView = imref2d(size(fixRefAdj));
end

% --- Apply transform to the entire moving cube ---
[rowsFix, colsFix] = size(fixRefAdj);
[~, ~, nBandsMov] = size(movCrop);
registeredCube = zeros(rowsFix, colsFix, nBandsMov, 'like', movCrop);

for b = 1:nBandsMov
    registeredCube(:,:,b) = imwarp(movCrop(:,:,b), tform, 'OutputView', outputView);
end

% --- Final overlay check ---
figure;
imshowpair(fixRefAdj, registeredCube(:,:,refBandMS), 'blend');
title('Overlay Fixed vs. Registered Moving — Band 9','FontSize',14);

disp('Registration complete. The MSI cube has been registered to the HSI reference.');
