function plotSupplyDemand(p_d,p_g,titlestring)
    
    % Determine the market clearing price which is where p_d and p_g intersect 
    % sort p_d in descending order
    [sorted_p_d, ~] = sort(p_d, 'descend');
    cumulative_demand = cumsum(ones(size(sorted_p_d)));
    
    % Sort p_g in ascending order
    [sorted_p_g, ~] = sort(p_g, 'ascend');
    cumulative_supply = cumsum(ones(size(sorted_p_g)));
    
    % Plotting the stairs plot for supply
    stairs(cumulative_supply, sorted_p_g, 'LineWidth', 2);
    hold on;
    
    % Plotting the stairs plot for demand
    stairs(cumulative_demand, sorted_p_d, 'LineWidth', 2);
    
    xlabel('Quantity');
    ylabel('Price');
    title(sprintf('%s\n',titlestring));
    legend('Supply', 'Demand');
    grid on
    
    % find the intersection of the two curves
    % find the index of the first element in sorted_p_d that is greater than sorted_p_g just look in the first 10 elements
    % this is the market clearing price
    for i = 1:15
        if sorted_p_d(i) < sorted_p_g(i)
            market_clearing_price = sorted_p_d(i);
            break;
        end
    end
    
    disp('Market clearing price: ');
    disp(market_clearing_price);

end