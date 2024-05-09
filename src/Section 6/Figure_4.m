theta = 0.2d0;

epsilon_values = 0.01:0.01:0.10;
num_epsilon_values = length(epsilon_values);

n_values = [2,3,4,5];
num_n_values = length(n_values);

num_trials = 1e6;
mu = -10;
line_styles = {'-', '--', ':', '-.'}; 
figure; % Create a new figure

for n_idx = 1:num_n_values
    n = n_values(n_idx);
    probability_return_values = zeros(1, num_epsilon_values);
    
    for epsilon_idx = 1:num_epsilon_values
        epsilon = epsilon_values(epsilon_idx);
        positive_diff_count = 0;
        true_diff_count = 0;

        rng(epsilon_idx * 10 + n_idx);  % Set the random seed based on epsilon and n indices

        for trial = 1:num_trials
            random_vector = mu + randn(1, n);

            squared_sum_diff = sum(random_vector)^2 - n * sum(random_vector.^2) + n * theta * theta;

            if squared_sum_diff >= 0
                positive_diff_count = positive_diff_count + 1;

                new_check = sum(random_vector) + sqrt(squared_sum_diff) - n * mu - sqrt(n) * (theta + norminv(1 - epsilon, 0, 1));
                if new_check >= 0
                    true_diff_count = true_diff_count + 1;
                end
            end
        end

        probability_return = true_diff_count / positive_diff_count;
        probability_return_values(epsilon_idx) = probability_return;
    end
    plot(epsilon_values, probability_return_values, 'LineWidth', 3, 'LineStyle', line_styles{n_idx}, 'DisplayName', ['n = ' num2str(n)]);
   
    hold on;
end

xlabel('$\varepsilon$', 'Interpreter', 'latex');
xlim([0.01 0.1])
% ylim([0.0 0.1])
ylabel('Fixing Probability');
% title('Probability Return vs. Epsilon for Different n');
grid on;
lgd = legend('Location', 'northeast', 'Interpreter', 'latex');
set(lgd, 'FontSize', 14);

