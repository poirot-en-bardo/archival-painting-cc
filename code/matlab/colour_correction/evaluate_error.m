%% Delta E Evaluation Function
function evaluate_error(ref_lab, corrected_lab, test_idx, m, n, space, output_folder, file_name)
    % Compute Delta E2000 Error
    deltaE2000_errors = deltaE2000(corrected_lab, ref_lab);
    
    % Calculate mean and max Delta E2000 for the test set
    mean_deltaE = mean(deltaE2000_errors(test_idx));
    max_deltaE = max(deltaE2000_errors(test_idx));
    
    % Compute the Delta E2000 error map for the full set (error map for entire dataset)
    error_map = reshape(deltaE2000(corrected_lab, ref_lab), m, n);
    
    % Display the results
    disp([space ' - Mean ΔE2000 Error: ', num2str(mean_deltaE)]);
    disp([space ' - Max ΔE2000 Error: ', num2str(max_deltaE)]);
    
    % Plot the error map
    figure;
    imagesc(error_map);
    colormap(jet);
    colorbar;
    clim([0 10]);  % Set the color axis limits for clarity
    title([space ' ΔE2000', ' (Mean: ', num2str(mean_deltaE, '%.2f'), ', Max: ', ...
        num2str(max_deltaE, '%.2f'), ')'], Interpreter="none");
    grid off;

    [rows, cols] = ind2sub([m, n], test_idx); % Convert linear indices to (row, col)
    hold on;  % Keep the error map displayed
    plot(cols, rows, 'w.', 'MarkerSize', 10);  % Plot white dots at training locations
    hold off;
    
    % Save the error map figure
    % saveas(gcf, fullfile(output_folder, space + file_name));
end