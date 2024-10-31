function [p_value, Hobs, df] = kruskalwallis_test(data)
    format long
    % Data cleaning
    % Identify and remove rows with any missing values
    na_values = any(ismissing(data), 2);
    data = data(~na_values, :); 

    % Initialize variables to store the indices of numeric and factor columns
    numeric_col = [];
    factor_col = [];
    
    % Loop to check  if the column is numeric or a factor 
    for col = 1:width(data)
        if isnumeric(data{:, col})
            numeric_col = col;
        elseif iscategorical(data{:, col}) || isstring(data{:, col}) || iscellstr(data{:, col})
            factor_col = col;
        end
    end
    
    % Check if numeric and factor columns were identified
    if isempty(numeric_col) || isempty(factor_col)
        error('Data must contain both a numeric column and a factor column.');
    end

    % Extract the numeric part to x and the factor part to g
    x = data{:, numeric_col};  
    g = data{:, factor_col};  
    
    % Convert grouping variable to cell array of strings if necessary
    if iscategorical(g)
        g = cellstr(g);
    elseif isstring(g)
        g = cellstr(g);
    end
    
    % Error checks
    if length(x) ~= length(g)
        error('Length of data vector and grouping variable must be the same.');
    end
    if length(unique(g)) < 3
        error('There must be at least three distinct groups in the grouping variable.');
    end
    
    % Calculation of ranks 
    ranks = tiedrank(x);  % This function calculates the rank of x with ties
    
    % Reorganization of ranks by groups
    unique_groups = unique(g);          % After calculating the rank of x these functions allocate 
    num_groups = length(unique_groups); % the rank of each value to its corresponding group
    ranks_cell = cell(1, num_groups);
    for i = 1:num_groups
        ranks_cell{i} = ranks(strcmp(g, unique_groups{i}));
    end
    
    % Total number of observations
    N = length(x);
    
    % Calculation of the tied values
    [~, ~, x_groups] = unique(x);
    tied = accumarray(x_groups, 1);
    
    % Correction factor z
    z = 1 - sum((tied.^3 - tied) / (N^3 - N)); % z is used to adjust the Hobs
    
    % Sum of ranks for each group
    sum_ranks = cellfun(@sum, ranks_cell);  % The sum of rank for each group has been calculated
    
    % Calculation of sum of squares of ranks divided by the number of observations in each group
    sum_ranks_squared_divided_by_n = sum((sum_ranks.^2) ./ cellfun(@length, ranks_cell));
    
    % Calculation of Hobs statistic
    Hobs = (12 / (N * (N + 1))) * sum_ranks_squared_divided_by_n - 3 * (N + 1);
    
    % Apply correction for ties to Hobs
    if any(tied > 1)
        Hobs = Hobs / z;
    end
    
    % Degrees of freedom
    df = num_groups - 1;
    
    % Calculation of p_value
    p_value = 1 - chi2cdf(Hobs, df);

    %Displaying results
    fprintf('p_value = %.10f\n df = %d\n H = %.2f\n', p_value, df, Hobs);
    if p_value > 0.05
        fprintf('H0 is accepted\n');
    else
        fprintf('H0 is rejected\n');
    end
end
