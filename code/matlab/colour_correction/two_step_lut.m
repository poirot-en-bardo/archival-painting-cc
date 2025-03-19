function two_step_color_correction_3dLUT()
    clc; clear; close all;
    roof = double(intmax('uint16')); % Normalization factor
    rng(10);

    %% ---------------- Load Hyperspectral Data ----------------
    cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_fuji-frame4.hdr";
    refFile  = "../../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

    hcube = hypercube(cubeFile);
    inCUBE = hcube.DataCube;
    bands = hcube.Wavelength;
    [m, n, bd] = size(inCUBE);

    hcube_ref = hypercube(refFile);
    refCUBE = hcube_ref.DataCube;
    bands_ref = hcube_ref.Wavelength;
    [mr, nr, bd_ref] = size(refCUBE);

    % Flatten the spectral data
    lincube = reshape(inCUBE, [], bd);
    lincube_ref = reshape(refCUBE, [], bd_ref);

    %% ---------------- Compute XYZ from Spectral Data ----------------
    % Load illuminant and CMFs
    ill = importdata("../../../data/CIE_D65.txt");
    CMFs_1931 = importdata("../../../data/CIE2degCMFs_1931.txt");
    CMFs = CMFs_1931;

    % For input cube:
    illIP = interp1(ill(:,1), ill(:,2), bands, 'spline');
    CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,3), bands, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,4), bands, 'spline')];
    sp_tristREF = CMFsIP .* illIP;
    xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);

    % For reference cube:
    illIP_ref = interp1(ill(:,1), ill(:,2), bands_ref, 'spline');
    CMFsIP_ref = [interp1(CMFs(:,1), CMFs(:,2), bands_ref, 'spline'), ...
                  interp1(CMFs(:,1), CMFs(:,3), bands_ref, 'spline'), ...
                  interp1(CMFs(:,1), CMFs(:,4), bands_ref, 'spline')];
    sp_tristREF_ref = CMFsIP_ref .* illIP_ref;
    xyz_ref = (lincube_ref * sp_tristREF_ref) ./ sum(sp_tristREF_ref(:,2), 1);

    %% ---------------- Reshape to 10x14 for Patch Indexing ----------------
    % Here we assume the chart is exactly 10 rows by 14 columns
    xyz_input_2d = reshape(xyz_input, 10, 14, 3);
    xyz_ref_2d   = reshape(xyz_ref,   10, 14, 3);

    %% ---------------- Step 1: LUT/Tone Correction Using Grayscale Patches ----------------
    % Define patches: border patches (first & last row/column) plus center (rows 5–6, cols 5–10)
    patchXYZ_in = [];
    patchXYZ_ref = [];

    % Border patches
    for row = 1:10
        for col = 1:14
            if row == 1 || row == 10 || col == 1 || col == 14
                patchVal_in = squeeze(xyz_input_2d(row, col, :))';
                patchVal_ref = squeeze(xyz_ref_2d(row, col, :))';
                patchXYZ_in  = [patchXYZ_in; patchVal_in];
                patchXYZ_ref = [patchXYZ_ref; patchVal_ref];
            end
        end
    end
    % Center region patches (rows 5–6, cols 5–10)
    for row = 5:6
        for col = 5:10
            patchVal_in = squeeze(xyz_input_2d(row, col, :))';
            patchVal_ref = squeeze(xyz_ref_2d(row, col, :))';
            patchXYZ_in  = [patchXYZ_in; patchVal_in];
            patchXYZ_ref = [patchXYZ_ref; patchVal_ref];
        end
    end

    % Build LUT correction for each channel using 1D interpolation
    correctedXYZ_2d_step1 = xyz_input_2d;
    for ch = 1:3
        inVals  = patchXYZ_in(:, ch);
        refVals = patchXYZ_ref(:, ch);
        [inVals_s, idx] = sort(inVals);
        refVals_s = refVals(idx);
        lutFun = @(x) interp1(inVals_s, refVals_s, x, 'linear', 'extrap');
        correctedXYZ_2d_step1(:,:,ch) = lutFun(xyz_input_2d(:,:,ch));
    end
    correctedXYZ_step1 = reshape(correctedXYZ_2d_step1, [], 3);

    %% ---------------- Convert Step 1 Corrected XYZ to RGB ----------------
    % Use built-in conversion; here we assume prophoto-rgb for consistency.
    rgb_corr1 = xyz2rgb(correctedXYZ_step1, 'ColorSpace', 'prophoto-rgb');
    rgb_ref   = xyz2rgb(xyz_ref, 'ColorSpace', 'prophoto-rgb');

    %% ---------------- Step 2: 3D LUT Correction using scatteredInterpolant ----------------
    % Use scatteredInterpolant to create a mapping from intermediate corrected RGB to reference RGB
    % (Both rgb_corr1 and rgb_ref are N-by-3 arrays; N = 10*14.)
    F_R = scatteredInterpolant(rgb_corr1(:,1), rgb_corr1(:,2), rgb_corr1(:,3), rgb_ref(:,1), 'linear', 'nearest');
    F_G = scatteredInterpolant(rgb_corr1(:,1), rgb_corr1(:,2), rgb_corr1(:,3), rgb_ref(:,2), 'linear', 'nearest');
    F_B = scatteredInterpolant(rgb_corr1(:,1), rgb_corr1(:,2), rgb_corr1(:,3), rgb_ref(:,3), 'linear', 'nearest');

    rgb_corr2 = zeros(size(rgb_corr1));
    for i = 1:size(rgb_corr1,1)
        rgb_corr2(i,1) = F_R(rgb_corr1(i,1), rgb_corr1(i,2), rgb_corr1(i,3));
        rgb_corr2(i,2) = F_G(rgb_corr1(i,1), rgb_corr1(i,2), rgb_corr1(i,3));
        rgb_corr2(i,3) = F_B(rgb_corr1(i,1), rgb_corr1(i,2), rgb_corr1(i,3));
    end

    %% ---------------- Reconstruct Final Corrected XYZ for Evaluation ----------------
    % Convert final corrected RGB (prophoto-rgb) back to XYZ
    correctedXYZ_step2 = rgb2xyz(rgb_corr2, 'ColorSpace', 'prophoto-rgb');

    %% ---------------- Evaluate & Display Results ----------------
    % For error evaluation, convert both reference and final corrected XYZ to Lab.
    lab_ref  = xyz2lab(reshape(xyz_ref, [], 3));
    lab_corr = xyz2lab(correctedXYZ_step2);
    
    % Here, we evaluate error over all pixels (or a train-test split if desired)
    test_idx = 1:size(lab_ref,1); % full dataset evaluation

    output_folder = '../../../results/error_maps';
    if ~exist(output_folder, 'dir'), mkdir(output_folder); end

    [~, img_name, ~] = fileparts(cubeFile);
    evaluate_error(lab_ref, lab_corr, test_idx, 10, 14, '2step_', output_folder, [img_name, '_2step_correction.png']);

    %% ---------------- Display Corrected Images ----------------
    % Convert original, step1, and step2 corrected XYZ to RGB for display
    rgb_original = xyz2rgb(reshape(xyz_input, 10, 14, 3), 'ColorSpace', 'prophoto-rgb'); % not ideal, but as a placeholder
    rgb_step1 = xyz2rgb(correctedXYZ_step1, 'ColorSpace', 'prophoto-rgb');
    rgb_step2 = xyz2rgb(correctedXYZ_step2, 'ColorSpace', 'prophoto-rgb');

    figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    nexttile;
    imagesc(reshape(rgb_original, 10, 14, 3));
    axis equal off;
    title('Original (RGB)');
    
    nexttile;
    imagesc(reshape(rgb_step1, 10, 14, 3));
    axis equal off;
    title('After Step1 LUT');
    
    nexttile;
    imagesc(reshape(rgb_step2, 10, 14, 3));
    axis equal off;
    title('After Step2 3D LUT');
    
    sgtitle('Two-Step Color Correction: LUT then 3D LUT');
end

%% -------------------------------------------------------------------------\n%% Polynomial Feature Function (unused in this version)\nfunction X_poly = poly3_features(input_data)\n    a = input_data(:,1);\n    b = input_data(:,2);\n    c = input_data(:,3);\n    feat1 = a; feat2 = b; feat3 = c;\n    feat4 = a.^2; feat5 = b.^2; feat6 = c.^2;\n    feat7 = a.*b; feat8 = a.*c; feat9 = b.*c;\n    feat10 = a.^3; feat11 = b.^3; feat12 = c.^3;\n    feat13 = a.^2.*b; feat14 = a.^2.*c; feat15 = b.^2.*a;\n    feat16 = b.^2.*c; feat17 = c.^2.*a; feat18 = c.^2.*b; feat19 = a.*b.*c;\n    X_poly = [feat1, feat2, feat3, feat4, feat5, feat6, feat7, feat8, feat9, ...\n              feat10, feat11, feat12, feat13, feat14, feat15, feat16, feat17, feat18, feat19];\nend\n"}
