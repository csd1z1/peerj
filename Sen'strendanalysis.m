% Read the base data to obtain projection and size information
[a, R] = geotiffread('F:\Data folder\2000.tif');
info = geotiffinfo('F:\Data folder\2000.tif');
[m, n] = size(a);

% Set the time span (2000-2022, 23 years in total)
cd = 2022 - 2000 + 1;  % Modify according to the actual year range
datasum = zeros(m*n, cd) + NaN;  % Initialize the data storage matrix

% Loop through and read multi-year data
k = 1;
for year = 2000:2022  % From start year to end year
    filename = ['F:\Data folder\', int2str(year), '.tif'];
    data = importdata(filename);
    data = reshape(data, m*n, 1);  % Convert to column vector
    datasum(:, k) = data;
    k = k + 1;
end

% Process invalid values (set values < -100 to NaN)
datasum(datasum < -100) = NaN;

% Initialize the result matrix
result = zeros(m, n) + NaN;

% Calculate Sen's slope for each pixel
for i = 1:size(datasum, 1)
    data = datasum(i, :);
    % Process only valid values (pixels with all values > 0)
    if min(data) > 0
        valuesum = [];
        % Calculate the rate of change for all year pairs
        for k1 = 2:cd
            for k2 = 1:(k1-1)
                cz = data(k1) - data(k2);  % Value change
                jl = k1 - k2;  % Time interval (years)
                value = cz ./ jl;  % Annual rate of change
                valuesum = [valuesum; value];
            end
        end
        % Use median as the Sen's slope result
        result(i) = median(valuesum);
    end
end

% Save the result as a GeoTIFF file
filename = ['F:\Data folder\Sentrendanalysis .tif'];
geotiffwrite(filename, result, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
disp('OK!');  % Indicate completion