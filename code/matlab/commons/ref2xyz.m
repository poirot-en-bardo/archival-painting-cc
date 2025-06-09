function XYZ = ref2xyz(S, cmf, refl)
    % S: [numBands x 1]         (illuminant)
    % cmf: [numBands x 3]       (color matching functions)
    % refl: [numColors x numBands]   (spectra in rows!)
    % Returns: XYZ [numColors x 3]

    k = 100 / dot(cmf(:,2), S); % CIE normalization

    % Compute XYZ for all colors
    % Each row of refl is a color; need to operate across columns
    X = k * sum(refl .* (cmf(:,1)' .* S'), 2);
    Y = k * sum(refl .* (cmf(:,2)' .* S'), 2);
    Z = k * sum(refl .* (cmf(:,3)' .* S'), 2);

    XYZ = [X, Y, Z];
end
