function [] = randPoints(num_of_points, upper_limit, lower_limit)
    
    range_power = log(upper_limit / lower_limit) / log(10);

    true_loc = zeros(num_of_points, 3);
    true_loc(:, 1) = lower_limit * 10 .^ (rand(num_of_points , 1) * range_power);
    true_loc(:, 2) = (rand(num_of_points , 1))*(2*pi);
    true_loc(:, 3) = (rand(num_of_points , 1))*(200);
    [true_loc(:, 1), true_loc(:, 2), true_loc(:, 3)] = pol2cart(true_loc(:, 2), true_loc(:, 1), true_loc(:, 3));

    save(fullfile("data/locations",'randLoc.mat'),'true_loc')

end