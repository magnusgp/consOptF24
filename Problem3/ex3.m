
% Assuming each p_d and p_g corresponds to a quantity of 1 unit
% Sorting demand prices in descending order for the demand curve
[sorted_p_d, ~] = sort(p_d, 'descend');
cumulative_demand = cumsum(ones(size(sorted_p_d)));

% Sorting generation prices in ascending order for the supply curve
[sorted_p_g, ~] = sort(p_g, 'ascend');
cumulative_supply = cumsum(ones(size(sorted_p_g)));

% Plotting the stairs plot for supply
figure;
stairs(cumulative_supply, sorted_p_g, 'LineWidth', 2);
hold on;

% Plotting the stairs plot for demand
stairs(cumulative_demand, sorted_p_d, 'LineWidth', 2);

% Adding labels and title
xlabel('Quantity');
ylabel('Price');
title('Supply-Demand Curve');
legend('Supply', 'Demand');

% find the market clearing price i.e. where supply is equal to demand




