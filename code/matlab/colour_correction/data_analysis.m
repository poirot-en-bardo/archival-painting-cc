
clc; clear; close all;
roof = double(intmax('uint16')); % Histogram normalization factor

% Load spectral data
ill = importdata('../../../data/CIE_D65.txt'); % Illuminant
CMFs_1931 = importdata('../../../data/CIE2degCMFs_1931.txt');
CMFs = CMFs_1931;
rng(10);

%% Select the reference and the cube to be corrected

cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_120f-velvia-f8.hdr";
% cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_fuji-frame4.hdr";
refFile = "../../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

%% Load and process the cubes
hcube = hypercube(cubeFile);
inCUBE = hcube.DataCube;
bands = hcube.Wavelength;
[m, n, bd] = size(inCUBE);
lincube = reshape(inCUBE, [], bd);

hcube_ref = hypercube(refFile);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;
[m, n, bd] = size(refCUBE);
lincube_ref = reshape(refCUBE, [], bd);

% Interpolate illuminant and CMFs
illIP = interp1(ill(:,1),ill(:,2),bands,'spline');
CMFsIP = [interp1(CMFs(:,1),CMFs(:,2),bands,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,3),bands,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,4),bands,'spline')];
sp_tristREF = CMFsIP .* illIP;

% Compute XYZ values for the input and reference cubes
xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2),1);

illIP = interp1(ill(:,1),ill(:,2),bands_ref,'spline');
CMFsIP = [interp1(CMFs(:,1),CMFs(:,2),bands_ref,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,3),bands_ref,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,4),bands_ref,'spline')];
sp_tristREF = CMFsIP .* illIP;

xyz_ref = (lincube_ref * sp_tristREF) ./ sum(sp_tristREF(:,2),1);

%% Correlation analysis
reshapedData = reshape(inCUBE, [], size(inCUBE, 3));
correlationMatrix = corrcoef(reshapedData);

imagesc(correlationMatrix);
colorbar;
title('Correlation Matrix of Spectral Bands');
xlabel('Spectral Bands');
ylabel('Spectral Bands');

%% PCA

[coeff, score, latent] = pca(reshapedData);

figure;
plot(mean(reshapedData(:, 1:30), 1)); % For the first block
hold on;
plot(mean(reshapedData(:, 31:60), 1)); % For the second block
plot(mean(reshapedData(:, 60:117), 1)); % For the third block
hold off;
legend({'Bands 1-30', 'Bands 31-60', 'Bands 60-117'});
title('Mean Spectral Response for Different Blocks');
xlabel('Wavelength (nm)');
ylabel('Mean Reflectance');


%%

Z = linkage(correlationMatrix, 'average');
dendrogram(Z);
title('Spectral Bands Hierarchical Clustering');
