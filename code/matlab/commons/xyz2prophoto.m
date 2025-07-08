function RGB = xyz2prophoto(XYZ, applyGamma)
    arguments
        XYZ (:,3) double
        applyGamma (1,1) logical = true
    end

    % 1) Linear ROMM-RGB
    M = [ 1.3460, -0.2556, -0.0511;
         -0.5446,  1.5082,  0.0205;
          0.0000,  0.0000,  1.2123 ];
    % simple N×3 → N×3 multiply
    RGB = XYZ * M.';

    % clamp negatives
    RGB = max(RGB, 0);

    % 2) Specular-gamma / toe if requested
    if applyGamma
        Et    = 1/512;
        slope = 16;
        gamma = 1.8;

        % apply piecewise
        mask = RGB < Et;
        RGB(mask)    = RGB(mask) * slope;
        RGB(~mask)   = RGB(~mask).^(1/gamma);
    end

    % 3) Clip to [0,1]
    RGB = min(RGB, 1);
end
