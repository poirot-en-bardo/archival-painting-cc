function XYZ = ref2xyz(S, cmf, relf)
    % S: [numBands x 1]
    % cmf: [numBands x 3]
    % d: [numBands x numPixels]
    % Returns: XYZ [numPixels x 3]
    k = 100 / dot(cmf(:,2), S); % Normalization using Y-bar
    % Each band for all pixels: [numBands x numPixels]
    X = k * sum(relf .* (cmf(:,1) .* S), 1);
    Y = k * sum(relf .* (cmf(:,2) .* S), 1);
    Z = k * sum(relf .* (cmf(:,3) .* S), 1);
    XYZ = [X; Y; Z]';
end