function XYZ_out = apply_msr_to_Y(XYZ_in, scales)
    % epsv = 1e-9;
    % 
    % % Extract luminance channel
    % Y = XYZ_in(:,:,2); 
    % 
    % % Normalize Y
    % Ynorm = max(Y ./ 100, epsv);
    % 
    % weight = 1 / numel(scales);
    % R = zeros(size(Ynorm));
    % 
    % % Multi-Scale Retinex on Y
    % for si = 1:numel(scales)
    %     s = scales(si);
    %     blur = imgaussfilt(Ynorm, s);
    %     R = R + weight * (log(Ynorm) - log(max(blur, epsv)));
    % end
    % 
    % % Rescale to [0,1]
    % R = R - min(R(:));
    % if max(R(:)) > 0
    %     R = R ./ max(R(:));
    % end
    % 
    % % Restore Y scale
    % Y_new = R * 100;
    % 
    % % Compute scaling factor and apply to all channels
    % scale = Y_new ./ (Y + epsv);
    % scale3 = repmat(scale, [1 1 3]);
    % XYZ_out = XYZ_in .* scale3;
    
end
