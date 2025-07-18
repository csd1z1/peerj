# peerj
# Time Series Data Analysis Toolbox  
This MATLAB toolbox implements three commonly used time series analysis methods: Gray Relational Analysis (GRA), Sen's Slope Estimator, and Mann-Kendall Trend Test. These methods help analyze and quantify long-term change trends in environmental, ecological, meteorological, and other domains.

## Functional Overview  
**1. Gray Relational Analysis (GRA)**  
-
-
-
-
-

**2. Sen's Slope Estimator**  
-
-
-
-

**3. Mann-Kendall Trend Test**  
-
-
-
-

## Installation & Usage  
**Environment Requirements**  
-
-

**Usage Procedure**  
1.
2.
3.

**Example Code**  
```
% Set base paths
dataPath = 'F:\Data folder\';
savePath = 'F:\Results folder\';

% Run Sen's Slope analysis
[a, R] = readgeoraster([dataPath, '2000.tif']);
info = geotiffinfo([dataPath, '2000.tif']);
[m, n] = size(a);

cd = 2022 - 2000 + 1;
datasum = zeros(m*n, cd) + NaN;

p = 1;
for year = 2000:2022
    filename = [dataPath, int2str(year), '.tif'];
    data = importdata(filename);
    data = reshape(data, m*n, 1);
    datasum(:, p) = data;
    p = p + 1;
end

result = zeros(m, n) + NaN;

for i = 1:size(datasum, 1)
    data = datasum(i, :);
    if min(data) > 0
        valuesum = [];
        for k1 = 2:cd
            for k2 = 1:(k1-1)
                cz = data(k1) - data(k2);
                jl = k1 - k2;
                value = cz ./ jl;
                valuesum = [valuesum; value];
            end
        end
        result(i) = median(valuesum);
    end
end

filename = [savePath, 'Sen_trend_result.tif'];
geotiffwrite(filename, result, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);

% Run Mann-Kendall trend test
sresult = zeros(m, n) + NaN;
for i = 1:size(datasum, 1)
    data = datasum(i, :);
    if min(data) > 0
        sgnsum = [];
        for k = 2:cd
            for j = 1:(k-1)
                sgn = sign(data(k) - data(j));
                sgnsum = [sgnsum; sgn];
            end
        end
        add = sum(sgnsum);
        sresult(i) = add;
    end
end

vars = cd * (cd-1) * (2*cd+5) / 18;
zc = zeros(m, n) + NaN;

sy = find(sresult == 0);
zc(sy) = 0;

sy = find(sresult > 0);
zc(sy) = (sresult(sy) - 1) ./ sqrt(vars);

sy = find(sresult < 0);
zc(sy) = (sresult(sy) + 1) ./ sqrt(vars);

geotiffwrite([savePath, 'MK_test_result.tif'], zc, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);

% Significance testing and result integration
data = importdata([savePath, 'MK_test_result.tif']);
sen_value = importdata([savePath, 'Sen_trend_result.tif']);
sen_value(abs(data) < 1.96) = NaN;

geotiffwrite([savePath, 'Significant_trends.tif'], sen_value, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);

% Run Gray Relational Analysis
clc;
clear;

file_path = [dataPath, 'your_data.xlsx'];
if exist(file_path, 'file')
    try
        data = readmatrix(file_path);
    catch
        data = xlsread(file_path);
    end
    
    [rows, cols] = size(data);
    
    if cols < 2
        error('Insufficient columns in data. At least 2 columns are required.');
    end
    
    data_normalized = mapminmax(data', 0.002, 1);
    data_normalized = data_normalized';
    
    figure(1)
    t = [1:rows];
    plot(t, data_normalized(:, 1), 'LineWidth', 2)
    hold on
    
    max_col = min(cols, 10);
    
    for i = 2:max_col
        plot(t, data_normalized(:, i), '--')
        hold on
    end
    
    grid on;
    xlabel('Samples');
    ylabel('Normalized Values');
    
    legend_entries = {'Y1'};
    for i = 2:max_col
        legend_entries{end+1} = ['x', num2str(i-1)];
    end
    legend(legend_entries, 'Location', 'best');
    title('Gray Relational Analysis - Normalized Data Comparison');
    
    gray_coefficients = zeros(1, max_col-1);
    disp('Gray relational coefficient calculation results:');
    
    for i = 2:max_col
        data_column = abs(data_normalized(:, i) - data_normalized(:, 1));
        d_max = max(data_column);
        d_min = min(data_column);
        a = 0.5;
        gray_relation = (d_min + a * d_max) ./ (data_column + a * d_max);
        coefficient = mean(gray_relation);
        gray_coefficients(i-1) = coefficient;
        disp(['x', num2str(i-1), ' and Y1 gray relational coefficient: ', num2str(coefficient)]);
    end
    
    figure(2)
    bar(gray_coefficients)
    grid on;
    xlabel('Variables');
    ylabel('Relational Coefficient');
    title('Comparison of Gray Relational Coefficients with Y1');
    
    x_ticks = 1:length(gray_coefficients);
    x_labels = cell(1, length(gray_coefficients));
    for i = 1:length(gray_coefficients)
        x_labels{i} = ['x', num2str(i)];
    end
    set(gca, 'XTick', x_ticks, 'XTickLabel', x_labels);
    
    [sorted_coeff, sorted_idx] = sort(gray_coefficients, 'descend');
    disp('Sorted by relational coefficient from highest to lowest:');
    for i = 1:length(sorted_coeff)
        disp(['x', num2str(sorted_idx(i)), ' coefficient: ', num2str(sorted_coeff(i))]);
    end
else
    error(['File does not exist: ', file_path]);
end
```

## Result Interpretation  
**Gray Relational Analysis**  
-
-
-

**Sen's Slope & Mann-Kendall Test**  
-
-
-

## Contribution & Feedback  
If you encounter issues or have improvement suggestions, please submit Issues or Pull Requests. We welcome all contributions!

## Citation  
If you use this toolbox in research, please appropriately cite the original methodology papers:

- **Gray Relational Analysis**  
Deng, J. L. (1982). Control problems of gray systems. Systems & Control Letters, 1(5), 288-294.

- **Sen's Slope Estimator**  
Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau. Journal of the American Statistical Association, 63(324), 1379-1389.

- **Mann-Kendall Test**  
Mann, H. B. (1945). Nonparametric tests against trend. Econometrica: Journal of the Econometric Society, 13(3), 245-259.

## License  
This project is licensed under the MIT License - see the LICENSE file for details.
