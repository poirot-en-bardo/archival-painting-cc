% --- Parameters ---
filename = '/Volumes/School/Thesis/Ageing/irradiance_plot.csv';   % your CSV file, with columns "x" and "y"
ref_wavelength = 340;    % nm
ref_irradiance  = 1.5;  % W/m^2/nm at the reference wavelength

% --- 1) Read the data ---
T = readtable(filename);    % requires headers "x" and "y"
x = T.x;
y = T.y;

% If your CSV has no header, you can instead do:
% data = readmatrix(filename);
% x = data(:,1);
% y = data(:,2);

% --- 2) Sort by increasing wavelength ---
[x_sorted, idx] = sort(x);
y_sorted        = y(idx);

% --- 3) Find y at the reference wavelength via interpolation ---
y_ref = interp1(x_sorted, y_sorted, ref_wavelength);

% --- 4) Compute scale factor so that E(ref_wavelength) = ref_irradiance ---
scale_factor = ref_irradiance / y_ref;

% --- 5) Apply scaling to get real spectral irradiance E(lambda) ---
E = y_sorted * scale_factor;

% --- 6) Integrate E(lambda) from min(x) to max(x) to get total irradiance ---
total_irradiance = trapz(x_sorted, E);

% --- 7) Display results ---
fprintf('Scale factor: %.4f\n', scale_factor);
fprintf('Total UV irradiance: %.2f W/m^2\n', total_irradiance);
