function RGB = xyz2prophoto(XYZ, applyGamma)
% Convert XYZ (D50) to ProPhoto RGB (ROMM RGB) with optional gamma correction
% XYZ: Nx3 matrix, rows are [X Y Z]
% applyGamma: set true to apply gamma 1.8 for display
% RGB: Nx3, [R G B] in ProPhoto RGB

    % ProPhoto RGB D50 transformation matrix (from Kodak/ICC)
    M = [  1.3459433, -0.2556075, -0.0511118;
          -0.5445989,  1.5081673,  0.0205351;
           0.0000000,  0.0000000,  1.2118128 ];

    if size(XYZ,2) ~= 3
        error('Input XYZ must be Nx3.');
    end

    % Linear ProPhoto RGB
    RGB = (M * XYZ')';

    % Clamp negatives to zero (optional, but recommended)
    RGB = max(0, RGB);

    % Apply gamma 1.8 if requested
    if nargin > 1 && applyGamma
        gamma = 1.8;
        RGB = RGB .^ (1/gamma);
    end

    % Clamp to [0,1] for display
    RGB = min(1, RGB);
end

