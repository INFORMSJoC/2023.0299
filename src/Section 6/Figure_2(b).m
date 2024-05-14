clear;
clc;

n_values = [2,3,4];
epsilon = 0.01:0.01:0.10; % Range of epsilon values

line_styles = {'-', '--', ':'}; 

figure;



for i = 1:length(n_values)
    n = n_values(i);
    y = 1-normcdf(sqrt(n)*norminv( (1 - epsilon).^(1/n),0,1 ),0,1);

    
    plot(epsilon, y, 'LineWidth', 3, 'LineStyle', line_styles{i}, 'DisplayName', ['n = ' num2str(n)]);
    hold on;
end

xlabel('$\varepsilon$', 'Interpreter', 'latex');
% ylabel('Value $1 - (1 - \varepsilon)^{1/n}$', 'Interpreter', 'latex');
xlim([0.01 0.1])
ylim([0.0 0.1])
ylabel('Fixing Probability');
% title('Portion of Violated Scenarios Identified by Corollary 4');
grid on;
lgd = legend('Location', 'Best', 'Interpreter', 'latex');

set(lgd, 'FontSize', 14);

