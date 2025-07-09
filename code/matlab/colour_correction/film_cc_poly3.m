clear; close all;

%% Load precomputed patch data
inputFile = '/home/oem/eliza/data/xyz_lab_rgb/colorchecker/yoda_halogen_fuji_underexp_balanced_colorchecker.mat';  
refFile = '/home/oem/eliza/masters-thesis/data/CCpassport_spectra.txt';
ill = importdata('../../../data/CIE_D50.txt');
CMFs = importdata('../../../data/CIE2degCMFs_1931.txt');

S = load(inputFile);

XYZ_input     = S.patchXYZ;
Lab_input     = S.patchLab;
RGB_input     = S.patchRGB;
RGB_lin_input = S.patchRGB_lin;


%%
% Load spectral reflectance reference file
refData = readmatrix(refFile, 'NumHeaderLines', 1);  
ref_wl = refData(:,1);   % Wavelengths (1st column)
ref_spectra = refData(:,2:end);    % 24 columns of reflectance data

valid_idx = ref_wl >= 380 & ref_wl <= 780;
ref_wl = ref_wl(valid_idx);
ref_spectra = ref_spectra(valid_idx, :);

% Interpolate illuminant and CMFs to match the wavelengths
illIP = interp1(ill(:,1), ill(:,2), ref_wl, 'spline');
CMFsIP = interp1(CMFs(:,1), CMFs(:,2:4), ref_wl, 'spline');  % x̄ ȳ z̄

% Combine to get the reference XYZ values
sp_tristREF = CMFsIP .* illIP;       % Element-wise multiply [N×3]
k = 100 / sum(sp_tristREF(:,2));     % Normalization so Y=100 for perfect white

patchXYZ_ref = k * (ref_spectra' * sp_tristREF);   % [24×N] × [N×3] = [24×3]
patchLab_ref      = xyz2lab_custom(patchXYZ_ref);
patchRGB_ref      = xyz2prophoto(patchXYZ_ref / 100, true);
patchRGB_lin_ref  = xyz2prophoto(patchXYZ_ref / 100, false);

XYZ_ref     = patchXYZ_ref;
Lab_ref     = patchLab_ref;
RGB_ref     = patchRGB_ref;
RGB_lin_ref = patchRGB_lin_ref;

%%

% Reshape to N×3
N = size(XYZ_input, 1);
m = 6; n = 4;  % for plotting

rng(42);  % for reproducibility
K = 6;
cv = cvpartition(N, 'KFold', K);

% Storage
mean_err = struct();
max_err  = struct();

fields = {'XYZ', 'Lab', 'RGB', 'RGB_lin'};
ref_struct = struct('XYZ', XYZ_ref, 'Lab', Lab_ref, 'RGB', RGB_ref, 'RGB_lin', RGB_lin_ref);
input_struct = struct('XYZ', XYZ_input, 'Lab', Lab_input, 'RGB', RGB_input, 'RGB_lin', RGB_lin_input);

% Loop over each representation
for f = 1:numel(fields)
    name = fields{f};
    fprintf('\nEvaluating correction in %s space\n', name);
    mean_de = zeros(K, 1);
    max_de  = zeros(K, 1);

    for k = 1:K
        train_idx = training(cv, k);
        test_idx  = test(cv, k);

        X_train = root_poly_features(input_struct.(name)(train_idx, :));
        % X_train = poly5_features(input_struct.(name)(train_idx, :));

        Y_train = ref_struct.(name)(train_idx, :);
        % mu_X = mean(X_train, 1);
        % sigma_X = std(X_train, 0, 1);
        % 
        % % Avoid division by zero
        % sigma_X(sigma_X == 0) = 1;
        
        % X_train = (X_train - mu_X) ./ sigma_X;

        % Linear regression
        coeffs = pinv(X_train) * Y_train;
        lambda = 0.001;
        % coeffs = (X_train' * X_train + lambda * eye(size(X_train, 2))) \ (X_train' * Y_train);
        coeffs = pinv(X_train' * X_train + lambda * eye(size(X_train,2))) * X_train' * Y_train;



        % Apply to all data
        X_all = root_poly_features(input_struct.(name));
        % X_all = poly5_features(input_struct.(name));

        Y_pred = X_all * coeffs;

        % X_all = (X_all - mu_X) ./ sigma_X;
        % Y_pred = X_all * coeffs;
        % Y_pred = Y_pred .* sigma_Y + mu_Y;  % inverse transform


        if strcmp(name, 'Lab')
            predLab = Y_pred;
            refLab = Lab_ref;
        
        elseif strcmp(name, 'RGB')
            % Convert predicted RGB to Lab via XYZ
            XYZ_pred = prophoto2xyz(Y_pred, true);
            predLab = xyz2lab_custom(XYZ_pred);
        
            refLab = Lab_ref;
        
        elseif strcmp(name, 'RGB_lin')
            XYZ_pred = prophoto2xyz(Y_pred, false);
            predLab = xyz2lab_custom(XYZ_pred);
        
            % Use spectral Lab ref
            refLab = Lab_ref;
        
        else  % XYZ
            predLab = xyz2lab_custom(Y_pred);
            refLab  = Lab_ref;
        end

        deltaE = deltaE2000(predLab, refLab);

        mean_de(k) = mean(deltaE(test_idx));
        max_de(k)  = max(deltaE(test_idx));
        % fprintf(' Fold %d | ΔE2000 Mean = %.2f | Max = %.2f\n', k, mean_de(k), max_de(k));
    end

    % Store results
    mean_err.(name) = mean_de;
    max_err.(name)  = max_de;
end

%%

% Normalize to [0, 1] for display
ref_RGB_display = patchRGB_ref;
ref_RGB_display = ref_RGB_display ./ max(ref_RGB_display(:));
input_RGB_display = RGB_input;
input_RGB_display = input_RGB_display ./ max(input_RGB_display(:));

figure;

subplot(1,2,1);
imshow(reshape(ref_RGB_display, [6, 4, 3])); axis off;
title('Reference RGB (from spectra)');

subplot(1,2,2);
imshow(reshape(input_RGB_display, [6, 4, 3])); axis off;
title('Camera RGB Input');

%% Summary
fprintf('\n---- Cross-Validation Summary (ΔE2000) ----\n');
for f = 1:numel(fields)
    name = fields{f};
    fprintf('%-8s | Mean = %.2f ± %.2f | Max = %.2f ± %.2f\n', ...
        name, ...
        mean(mean_err.(name)), std(mean_err.(name)), ...
        mean(max_err.(name)),  std(max_err.(name)));
end

%% --- 3rd-degree polynomial features
function X_poly = poly3_features(input_data)
    a = input_data(:,1);
    b = input_data(:,2);
    c = input_data(:,3);
    X_poly = [ ...
        ones(size(a)), ...
        a, b, c, ...
        a.^2, b.^2, c.^2, ...
        a.*b, a.*c, b.*c, ...
        a.^3, b.^3, c.^3, ...
        a.^2.*b, a.^2.*c, ...
        b.^2.*a, b.^2.*c, ...
        c.^2.*a, c.^2.*b, ...
        a.*b.*c];
end

function X_rootpoly = root_poly_features(input_data)
    input_data = normalize(input_data, 'range');  % Scales to [0,1]


    a = input_data(:,1);  % R or channel 1
    b = input_data(:,2);  % G or channel 2
    c = input_data(:,3);  % B or channel 3

    % Avoid negative roots and division by zero
    a = max(a, 0); b = max(b, 0); c = max(c, 0);

    X_rootpoly = [ ...
        ones(size(a)), ...
        a, b, c, ...
        sqrt(a), sqrt(b), sqrt(c), ...
        a.*b, a.*c, b.*c, ...
        sqrt(a.*b), sqrt(a.*c), sqrt(b.*c), ...
        a.^2, b.^2, c.^2, ...
        sqrt(abs(a.^2 + b.^2 + c.^2))];
    % X_rootpoly = [ ...
    % ones(size(a)), ...
    % a, b, c, ...
    % sqrt(a), sqrt(b), sqrt(c), ...
    % a.*b, a.*c, b.*c, ...
    % sqrt(a.*b), sqrt(a.*c), sqrt(b.*c), ...
    % a.^2, b.^2, c.^2, ...
    % sqrt(abs(a.^2 + b.^2 + c.^2)), ...
    % a ./ (b + eps), b ./ (c + eps), c ./ (a + eps)];  % avoid div by zero

end

function X_poly = poly4_features(input_data)
    a = input_data(:,1);
    b = input_data(:,2);
    c = input_data(:,3);

    X_poly = [ ...
        ones(size(a)), ...                              % Degree 0
        a, b, c, ...                                     % Degree 1
        a.^2, b.^2, c.^2, a.*b, a.*c, b.*c, ...          % Degree 2
        a.^3, b.^3, c.^3, ...                            % Degree 3
        a.^2.*b, a.^2.*c, b.^2.*a, b.^2.*c, c.^2.*a, c.^2.*b, ...
        a.*b.*c, ...                                     % Mixed degree-3
        a.^4, b.^4, c.^4, ...                            % Degree 4
        a.^3.*b, a.^3.*c, b.^3.*a, b.^3.*c, c.^3.*a, c.^3.*b, ...
        a.^2.*b.^2, a.^2.*c.^2, b.^2.*c.^2, ...
        a.^2.*b.*c, a.*b.^2.*c, a.*b.*c.^2 ...           % Full mixed degree-4 terms
    ];
end

function X_poly = poly5_features(input_data)
    a = input_data(:,1);
    b = input_data(:,2);
    c = input_data(:,3);

    X_poly = [ ...
        ones(size(a)), ...                             % Degree 0
        a, b, c, ...                                    % Degree 1
        a.^2, b.^2, c.^2, a.*b, a.*c, b.*c, ...         % Degree 2
        a.^3, b.^3, c.^3, ...
        a.^2.*b, a.^2.*c, b.^2.*a, b.^2.*c, c.^2.*a, c.^2.*b, ...
        a.*b.*c, ...                                    % Degree 3 mixed
        a.^4, b.^4, c.^4, ...
        a.^3.*b, a.^3.*c, b.^3.*a, b.^3.*c, c.^3.*a, c.^3.*b, ...
        a.^2.*b.^2, a.^2.*c.^2, b.^2.*c.^2, ...
        a.^2.*b.*c, a.*b.^2.*c, a.*b.*c.^2, ...         % Degree 4 mixed
        a.^5, b.^5, c.^5, ...
        a.^4.*b, a.^4.*c, b.^4.*a, b.^4.*c, c.^4.*a, c.^4.*b, ...
        a.^3.*b.^2, a.^3.*c.^2, b.^3.*a.^2, b.^3.*c.^2, c.^3.*a.^2, c.^3.*b.^2, ...
        a.^3.*b.*c, a.*b.^3.*c, a.*b.*c.^3, ...
        a.^2.*b.^2.*c, a.^2.*b.*c.^2, a.*b.^2.*c.^2 ...  % Full mixed degree-5 terms
    ];
end


function X_poly = poly2_features(input_data)
    a = input_data(:,1);
    b = input_data(:,2);
    c = input_data(:,3);

    X_poly = [ ...
        ones(size(a)), ...          % Constant term (bias)
        a, b, c, ...                % Linear terms
        a.^2, b.^2, c.^2, ...       % Quadratic terms
        a.*b, a.*c, b.*c ...        % Interaction terms
    ];
end

function X_lin = linear_features(input_data)
    a = input_data(:,1);
    b = input_data(:,2);
    c = input_data(:,3);

    X_lin = [ ...
        ones(size(a)), ...  % Bias term (optional, but standard in linear regression)
        a, b, c ...         % Linear terms
    ];
end

