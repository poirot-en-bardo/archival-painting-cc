function XYZ = prophoto2xyz(RGB, applyGamma)
    arguments
        RGB {mustBeNumeric}
        applyGamma (1,1) logical = true
    end

    original_shape = size(RGB);

    % Reshape if input is an image (HxWx3)
    if ndims(RGB) == 3 && original_shape(3) == 3
        RGB = reshape(RGB, [], 3);
        reshape_back = true;
    elseif size(RGB,2) ~= 3
        error('Input must have 3 channels (Nx3 or HxWx3).');
    else
        reshape_back = false;
    end

    % 1) Inverse gamma correction
    if applyGamma
        Et    = 1/512;   % Match xyz2prophoto
        slope = 16;
        gamma = 1.8;

        mask = RGB < Et * slope;
        RGB(mask)  = RGB(mask) / slope;
        RGB(~mask) = RGB(~mask) .^ gamma;
    end

    % 2) Linear ProPhoto to XYZ
    M_inv = [ ...
        0.7976749, 0.1351917, 0.0313534;
        0.2880402, 0.7118741, 0.0000857;
        0.0000000, 0.0000000, 0.8252100 ];
    
    XYZ = 100 .* (RGB * M_inv.');

    % Reshape back to image if needed
    if reshape_back
        XYZ = reshape(XYZ, original_shape);
    end
end
