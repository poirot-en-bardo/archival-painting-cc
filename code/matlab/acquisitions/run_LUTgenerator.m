%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script performs the color correction of an input image
% to be as close as possible to a reference image.
% The color corrected image is written on disk and the corresponding
% 3D-LUT is provided as output.
%
%                                       Written by Giorgio Trumpy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;

maxDeltaE = 16; % range of the residual deltaE image

%% Read, crop and register the images

f = msgbox('Select the image to be corrected');
movegui(f,'north')

[fileIN,pathIN] = uigetfile('*.tif');

close(f)

f = msgbox(sprintf(['Select the reference image for \n' fileIN]));
movegui(f,'north')

[fileREF,pathREF] = uigetfile([pathIN '*.tif']);

close(f)

roof = double(intmax('uint16'));

inputIMGfull = imread([pathIN fileIN]);
refIMG = imread([pathREF fileREF]);

answerCROP = questdlg('Do the images need to be copped?', ...
'PRE-CROP','Yes','No','Yes');

switch answerCROP
    case 'Yes'

    % unnecessary parts are excluded by cropping (e.g. perforation)
    imshow(inputIMGfull)
    title('drag-and-drop inside the image part');
    h = imrect;
    posIP = getPosition(h);
    close

    % unnecessary parts are excluded by cropping (e.g. perforation)
    imshow(refIMG)
    title('drag-and-drop inside the image part');
    h = imrect;
    posRF = getPosition(h);
    close

    inputIMGcrop = inputIMGfull(round(posIP(2)):round(posIP(2)+posIP(4)), ...
        round(posIP(1)):round(posIP(1)+posIP(3)),:);
    refIMGcrop = refIMG(round(posRF(2)):round(posRF(2)+posRF(4)), ...
        round(posRF(1)):round(posRF(1)+posRF(3)),:);

    case 'No'

    inputIMGcrop = inputIMGfull;
    refIMGcrop = refIMG;

end

% Registration
answerREG = questdlg('Do the images need to be registered?', ...
'REGISTRATION','Yes','No','Yes');

switch answerREG
    case 'Yes'

    Moving_Image = rgb2gray(inputIMGcrop);
    Fixed_Image = rgb2gray(refIMGcrop);
    registrationEstimator % launch the Registration GUI
    f = warndlg(['Load images from workspace. Export registered image ' ...
        'with the default name. When finished press OK.'], 'Registration');
    waitfor(f);

    % The transformation found is applied to the input image
    inputIMGreg = uint16(zeros(size(Fixed_Image,1),size(Fixed_Image,2),3));
    fixedRef = imref2d(movingReg.SpatialRefObj.ImageSize);
    for i=1:3
        inputIMGreg(:,:,i) = imwarp(inputIMGcrop(:,:,i),...
            movingReg.Transformation,'OutputView',fixedRef);
    end

    % The black pixel at the edge of the registered image are cropped out
    g = rgb2gray(inputIMGreg);
    g = [g(1,:) g(end,:) g(:,1)' g(:,end)'];
    while any(g==0)
        inputIMGreg = inputIMGreg(2:(end-1),2:(end-1),:);
        g = rgb2gray(inputIMGreg);
        g = [g(1,:) g(end,:) g(:,1)' g(:,end)'];
    end
    crop = (size(refIMGcrop,1)-size(inputIMGreg,1))/2;
    refIMGcrop = refIMGcrop(crop+1:end-crop,crop+1:end-crop,:);

    case 'No'
        
    if strcmp(answerCROP, 'No')       
        inputIMGreg = inputIMGcrop;
    else
        inputIMGreg = ...
            inputIMGfull(round(posRF(2)):round(posRF(2)+posRF(4)), ...
            round(posRF(1)):round(posRF(1)+posRF(3)),:);
    end

end

clear movingReg

%% Find the LUTs for the transformation

[inputSEG,refSEG,LUT_R,LUT_G,LUT_B] = findLUT(inputIMGreg,refIMGcrop);

%% Display result

disp('Displaying results...')

INval = (0:1:roof)';
figure,plot(INval,LUT_R,'r','LineWidth',2)
axis([0 roof 0 roof])
hold on
plot(INval,LUT_G,'g','LineWidth',2)
plot(INval,LUT_B,'b','LineWidth',2)

luttedIMGreg = applyLUT(inputIMGreg,LUT_R,LUT_G,LUT_B);
figure,imshowpair(luttedIMGreg,refIMGcrop,'checkerboard');

dE2Kimg = imcolordiff(refIMGcrop,luttedIMGreg,'Standard',"CIEDE2000");
figure,imshow(dE2Kimg,'DisplayRange',[0 maxDeltaE],'Colormap',jet(255))
pause(0.2)
colorbar
title('residual deltaE');
pause(0.2)

luttedIMGfull = applyLUT(inputIMGfull,LUT_R,LUT_G,LUT_B);
figure,imshow(luttedIMGfull);

imwrite(luttedIMGfull,[pathIN erase(fileIN,'.tif') '_wLUT.tif'], 'tif')

%% Save 1D-LUT in standard format

rgbLUT = [LUT_R LUT_G LUT_B];

oneDlut = [round(rgbLUT(4:4:end,:)/4) (1:16384)'];

dlmwrite([pathIN erase(fileIN,'.tif') '_pos2ref_1D-LUT.ilut'],oneDlut);

%% end

clc
