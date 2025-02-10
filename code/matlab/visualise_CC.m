%% Define some parameters
clc; clear; close all;
roof = double(intmax('uint16'));

histoFACT = 200; % Factor used to find the edges of the histograms

ill = importdata('../data/CIE_D65.txt');
CMFs = importdata('../data/CIE2degCMFs.txt');

%% Select input folder containing HDR cubes
f = msgbox('Select the folder containing HDR cubes');
movegui(f, 'north');
pathCB = uigetdir();
close(f);

% Get a list of all .hdr files in the folder
hdrFiles = dir(fullfile(pathCB, '*.hdr'));

% Create output folder if it doesn't exist
outputFolder = fullfile(pathCB, 'output_images');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Process each HDR file
for i = 1:length(hdrFiles)
    hdrFileName = hdrFiles(i).name;
    fullHdrPath = fullfile(pathCB, hdrFileName);

    % Load hyperspectral cube
    hcube = hypercube(fullHdrPath);
    inCUBE = hcube.DataCube;
    bands = hcube.Wavelength;
    [m, n, bd] = size(inCUBE);
    lincube = reshape(inCUBE, [], bd);

    % Interpolate illuminant and CMFs to captured wavelengths
    illIP = interp1(ill(:,1), ill(:,2), bands, 'spline');
    CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,3), bands, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,4), bands, 'spline')];
    sp_tristREF = CMFsIP .* illIP;

    % Calculate XYZ tristimulus values
    trist = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2),1);
    rgb = xyz2rgb(trist, 'ColorSpace', 'adobe-rgb-1998');
    inputIMG = uint16(reshape(rgb, m, n, 3) .* roof);

    % Extract title from HDR filename (everything after the first "_" before the extension)
    underscoreIdx = strfind(hdrFileName, '_');
    if isempty(underscoreIdx)
        titleText = hdrFileName(1:end-4); % Remove ".hdr" if no underscore
    else
        titleText = hdrFileName(underscoreIdx(1)+1:end-4); % Remove ".hdr"
    end
    titleText = strrep(titleText, '_', ' '); % Replace underscores with spaces in title

    %% Create an image with a title overlay
    % Convert input image to double for annotation
    overlayIMG = im2double(imresize(inputIMG, 100, 'nearest'));

    % Create a figure and add the title
    figure('Visible', 'off'); % Hide figure window for processing
    imshow(overlayIMG);
    hold on;
    text(10, 30, ['Visualisation: ', titleText], 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
    hold off;

    % Save the image with the title overlay
    outputImagePath = fullfile(outputFolder, [titleText, '.png']);
    frame = getframe(gca);
    imwrite(frame.cdata, outputImagePath);

    % Close the figure after saving
    close;
end

disp('Processing complete. PNG images saved in output_images folder.');
