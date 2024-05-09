% Clear the command window and workspace
clear;
clc;

% Define values of n and epsilon
n_values = [2, 3, 4];
eps = 0.1;

% Define line styles for different n values
line_styles = {'-', '--', ':'}; 

% Create a new figure
figure;

% Loop through each value of n
for i = 1:length(n_values)
    n = n_values(i);
    
    % Calculate the range of p values
    p_range = 1 - (eps^(1/n)):0.01:1;
    
    % Calculate y values using the updated formula
    y = (1 - p_range).^n;

    % Plot the data for the current n value
    plot(p_range, y, 'LineWidth', 3, 'LineStyle', line_styles{i}, 'DisplayName', ['n = ' num2str(n)]);
    
    % Hold the plot to add more data
    hold on;
end

% Label the axes and add title
xlabel('$\hat{p}$', 'Interpreter', 'latex');
% ylabel('Value $(1 - p)^n$', 'Interpreter', 'latex');
% title('Value $(1 - p)^n$ vs. $p$');
ylim([0.0 0.1])
ylabel('Fixing Probability');
% Display grid lines and legend
grid on;
% legend;
lgd = legend('Location', 'Best', 'Interpreter', 'latex');

set(lgd, 'FontSize', 14);


