function RGB = xyz2prophoto(XYZ, applyGamma)
    arguments
        XYZ {mustBeNumeric}
        applyGamma (1,1) logical = true
    end

    original_shape = size(XYZ);

    % Reshape if input is an image (HxWx3)
    if ndims(XYZ) == 3 && original_shape(3) == 3
        XYZ = reshape(XYZ, [], 3);
        reshape_back = true;
    elseif size(XYZ,2) ~= 3
        error('Input must have 3 channels (Nx3 or HxWx3).');
    else
        reshape_back = false;
    end

    % 1) Linear ROMM-RGB conversion
    M = [ ...
         1.3459433, -0.2556075, -0.0511118;
        -0.5445989,  1.5081673,  0.0205351;
         0.0000000,  0.0000000,  1.2118128 ];
    RGB = XYZ * M.';

    % Clamp negatives before gamma (for stability)
    RGB = max(RGB, 0);

    % 2) Gamma encoding
    if applyGamma
        Et    = 1/512;
        slope = 16;
        gamma = 1.8;

        mask = RGB < Et;
        RGB(mask)  = RGB(mask) * slope;
        RGB(~mask) = RGB(~mask).^(1/gamma);
    end

    % 3) Final clipping to [0, 1] for display
    RGB = min(RGB, 1);

    if reshape_back
        RGB = reshape(RGB, original_shape);
    end
end
