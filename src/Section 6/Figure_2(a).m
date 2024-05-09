clear;
clc;

epsilon_values = [0.01, 0.04, 0.07, 0.10];
n_values = 1:1:5; % Range of n values

line_styles = {'-', '--', ':', '-.'}; 

figure;

for i = 1:length(epsilon_values)
    epsilon = epsilon_values(i);
    y = 1 - normcdf(sqrt(n_values) .* norminv((1 - epsilon).^(1./n_values), 0, 1), 0, 1);

    plot(n_values, y, 'LineWidth', 3, 'LineStyle', line_styles{i}, 'DisplayName', ['$\varepsilon = ' num2str(epsilon, '%.2f') '$']);
    hold on;
end

xlabel('$n$', 'Interpreter', 'latex');
% ylabel('Value $1 - (1 - \varepsilon)^{1/n}$', 'Interpreter', 'latex');
% title('Portion of Violated Scenarios Identified by Corollary 4');
grid on;

xticks(n_values);  % Set custom x-axis tick values
xticklabels(n_values);  % Set custom x-axis tick labels

ylim([0.0 0.1])
ylabel('Fixing Probability');


% Place the legend in the best location and use LaTeX interpreter
lgd = legend('Location', 'Best', 'Interpreter', 'latex');

set(lgd, 'FontSize', 14);


