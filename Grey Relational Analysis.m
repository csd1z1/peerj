clc;
clear;

% Read data
file_path = '';
if exist(file_path, 'file')
    try
        % Compatible with newer MATLAB versions
        data = readmatrix(file_path);
    catch
        % Fallback for older MATLAB versions
        data = xlsread(file_path);
    end
    
    [rows, cols] = size(data);
    
    % Display basic data information
    disp(['Data dimensions: ', num2str(rows), ' rows x ', num2str(cols), ' columns']);
    
    % Check if data is sufficient
    if cols < 2
        error('Insufficient columns in data. At least 2 columns are required.');
    end
    
    % Data normalization
    data_normalized = mapminmax(data', 0.002, 1);
    data_normalized = data_normalized';
    
    % Plot normalized data comparison
    figure(1)
    t = [1:rows];
    plot(t, data_normalized(:, 1), 'LineWidth', 2)
    hold on
    
    % Dynamically adjust loop range to ensure it does not exceed column count
    max_col = min(cols, 10); % Display up to 10 columns
    
    for i = 2:max_col
        plot(t, data_normalized(:, i), '--')
        hold on
    end
    
    grid on;
    xlabel('Samples');
    ylabel('Normalized Values');
    
    % Dynamically create legend
    legend_entries = {'Y1'};
    for i = 2:max_col
        legend_entries{end+1} = ['x', num2str(i-1)];
    end
    legend(legend_entries, 'Location', 'best');
    title('Grey Relational Analysis - Normalized Data Comparison');
    
    % Calculate grey relational coefficients
    grey_coefficients = zeros(1, max_col-1);
    disp('Grey relational coefficient calculation results:');
    
    for i = 2:max_col
        % Calculate absolute difference from reference column
        data_column = abs(data_normalized(:, i) - data_normalized(:, 1));
        
        % Get global maximum and minimum values
        d_max = max(data_column);
        d_min = min(data_column);
        
        % Calculate grey relational matrix
        a = 0.5; % Resolution coefficient
        grey_relation = (d_min + a * d_max) ./ (data_column + a * d_max);
        coefficient = mean(grey_relation);
        grey_coefficients(i-1) = coefficient;
        
        disp(['x', num2str(i-1), ' and Y1 grey relational coefficient: ', num2str(coefficient)]);
    end
    
    % Plot relational coefficient bar chart
    figure(2)
    bar(grey_coefficients)
    grid on;
    xlabel('Variables');
    ylabel('Relational Coefficient');
    title('Comparison of Grey Relational Coefficients with Y1');
    
    % Add variable labels
    x_ticks = 1:length(grey_coefficients);
    x_labels = cell(1, length(grey_coefficients));
    for i = 1:length(grey_coefficients)
        x_labels{i} = ['x', num2str(i)];
    end
    set(gca, 'XTick', x_ticks, 'XTickLabel', x_labels);
    
    % Display sorted relational coefficients
    [sorted_coeff, sorted_idx] = sort(grey_coefficients, 'descend');
    disp('Sorted by relational coefficient from highest to lowest:');
    for i = 1:length(sorted_coeff)
        disp(['x', num2str(sorted_idx(i)), ' coefficient: ', num2str(sorted_coeff(i))]);
    end
else
    error(['File does not exist: ', file_path]);
end