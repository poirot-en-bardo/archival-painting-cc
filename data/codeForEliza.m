%% Define some parameters
clc
roof = double(intmax('uint16'));

histoFACT = 200; % factor used to find the edges of the histograms

ill = importdata('CIE_D65.txt');
CMFs = importdata('CIE2degCMFs.txt');

%% Select input cube and calculate input RGB image

f = msgbox('Select the cube to be corrected');
movegui(f,'north')
[headIN,pathCB] = uigetfile('*.hdr');
close(f)
hcube = hypercube([pathCB headIN]);
inCUBE = hcube.DataCube;
bands = hcube.Wavelength;
[m,n,bd] = size(inCUBE);

lincube = reshape(inCUBE,[],bd);

% Interpolate illuminant and CMF to the captured wavelengths
illIP = interp1(ill(:,1),ill(:,2),bands,'spline');
CMFsIP = [interp1(CMFs(:,1),CMFs(:,2),bands,'spline') ...
    interp1(CMFs(:,1),CMFs(:,3),bands,'spline') ...
    interp1(CMFs(:,1),CMFs(:,4),bands,'spline')];
sp_tristREF = CMFsIP.*illIP;

%calculate XYZ tristimulus values
trist = (lincube * sp_tristREF)./sum(sp_tristREF(:,2),1);
rgb = xyz2rgb(trist,'ColorSpace','adobe-rgb-1998');
inputIMG = uint16(reshape(rgb,m,n,3).*roof);

%% Display input image
%resize for display as otherwise the image is very small, only 140 pixels
figure, imshow(imresize(inputIMG, 100, 'nearest'));