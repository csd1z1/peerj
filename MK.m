% Set base paths (modify according to your actual situation)
dataPath = 'F:\Data folder\';
savePath = 'F:\Data folder\';

% Read base data to obtain projection and size information
[a, R] = readgeoraster([dataPath, '2000.tif']);
info = geotiffinfo([dataPath, '2000.tif']);
[m, n] = size(a);

% Set time span (2000-2022, 23 years in total)
cd = 2022 - 2000 + 1;
datasum = zeros(m*n, cd) + NaN;

% Loop through and read multi-year data
p = 1;
for year = 2000:2022
    filename = [dataPath, int2str(year), '.tif'];
    data = importdata(filename);
    data = reshape(data, m*n, 1);
    datasum(:, p) = data;
    p = p + 1;
end

% Perform Mann-Kendall trend test
sresult = zeros(m, n) + NaN;
for i = 1:size(datasum, 1)
    data = datasum(i, :);
    % Process only valid values (all values > 0)
    if min(data) > 0
        sgnsum = [];
        % Calculate change directions for all year pairs
        for k = 2:cd
            for j = 1:(k-1)
                sgn = sign(data(k) - data(j));  % Simplified sign function
                sgnsum = [sgnsum; sgn];
            end
        end
        % Sum all change directions
        add = sum(sgnsum);
        sresult(i) = add;
    end
end

% Calculate standard normal distribution statistic Zc
vars = cd * (cd-1) * (2*cd+5) / 18;
zc = zeros(m, n) + NaN;

% Compute Zc values based on sresult
sy = find(sresult == 0);
zc(sy) = 0;

sy = find(sresult > 0);
zc(sy) = (sresult(sy) - 1) ./ sqrt(vars);

sy = find(sresult < 0);
zc(sy) = (sresult(sy) + 1) ./ sqrt(vars);

% Save MK test results
geotiffwrite([savePath, 'MK_test_results.tif'], zc, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);

% Significance test and integrate with Sen's slope results
[a, R] = readgeoraster([dataPath, '2000.tif']);
info = geotiffinfo([dataPath, '2000.tif']);
data = importdata([savePath, 'MK_test_results.tif']);
sen_value = importdata([savePath, 'Sen_based_GPP_trend.tif']);

% Retain only pixels passing the 95% significance test
sen_value(abs(data) < 1.96) = NaN;

% Save final results
geotiffwrite([savePath, 'GPP_MK_Sen_significant_trends.tif'], sen_value, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);

disp('Process completed!');