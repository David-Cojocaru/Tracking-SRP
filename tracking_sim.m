
% Parameters
dt = 0.05;    % Sampling time
fps = 1 / dt;   % Frames per second
T = 2000;  % Total number of time steps
dim = 3;

% Initialization - here object 1 is the "seeker", while object 2 is the moving target
x1 = zeros(T, dim) + [150, 0, 0];
x2 = zeros(T, dim) + [0, 0, 0];
v1 = zeros(T, dim) + [0, 0, 0];
v2 = zeros(T, dim) + [0, -5, 0];

x_err = zeros(T, 1);
v_err = zeros(T, 1);
SNR_arr = zeros(T, 1);
dist_err = zeros(T, 1);
DOA_err = zeros(T, 1);

u = [0, 0, 0];
u_limit = 25; % limited to an acceleration of 50 m/s^2
v_limit = 50; % limited to a velocity of 50 m/s
%u2 = (2*rand(T, dim) - 1) * 5; % we generate random acceleration step for each time step to get a semi-random target route

dist_goal = 5; % the desired final relative distance
vel_goal = 5; % the desired final relative velocity

sum_error = 0;
prev_error = 0;

prepro = 1;
I = 0;
estim_dist_mode = 1;
estim_DOA_mode = 1;
r0 = 100;

x_dif_estim_dist = 0;
fig_ind = 3;

%%
K_p = 500 * 10^-3;  % Proportional gain
K_i = 20 * 10^-3;  % Integral gain
K_d = 1000 * 10^-3;  % Derivative gain

% Control loop
for i = 1:T
    
    % Apply control input u1 to the system (e.g., update acceleration of object 1)
    if i > 1
        
        v1(i, :) = v1(i-1, :) + u * dt;
        
        if norm(v1) > v_limit
            disp(['step ' num2str(i) ' over limit'])
            v1 = v1 * v_limit / norm(v1);
        end
        
        v2(i, :) = v2(i-1, :) + u2(i, :) * dt;
        x1(i, :) = x1(i-1, :) + v1(i, :) * dt;
        x2(i, :) = x2(i-1, :) + v2(i, :) * dt;
        
        prepro = 0;
    end
    
    % Measure relative position and velocity
    x_dif = x2(i, :) - x1(i, :);
    v_dif = v2(i, :) - v1(i, :);
    
    x_dif_true_dist = norm(x_dif);
    x_dif_true_DOA = x_dif / norm(x_dif);
    
    if estim_DOA_mode
        if prepro
            [SNR, x_dif_estim_DOA, Interp] = calc_DOA_SRP(x_dif, i, r0, prepro);
        else
            [SNR, x_dif_estim_DOA] = calc_DOA_SRP(x_dif, i, r0, prepro);
        end
    end
    
    if estim_dist_mode
        x_dif_estim_dist = Interp(SNR);
    end
    
    
    
    if estim_dist_mode
        dist_err(i) = abs(x_dif_estim_dist - norm(x_dif));
    end
    
    if estim_DOA_mode
        DOA_err(i) = rad2deg(acos(...
            (x_dif_estim_DOA*transpose(x_dif_true_DOA))./(sqrt(sum(x_dif_estim_DOA.^2, 2))*norm(x_dif_true_DOA))...
            ));
        
        SNR_arr(i) = SNR;
    end

    % Compute error
    if estim_dist_mode
        if estim_DOA_mode
            error = x_dif_estim_dist * x_dif_estim_DOA;  % Use position error or v_dif for velocity error
        else
            error = x_dif_estim_dist * x_dif_true_DOA;
        end
    else
        if estim_DOA_mode
            error = x_dif_true_dist * x_dif_estim_DOA;  % Use position error or v_dif for velocity error
        else
            error = x_dif;
        end
    end
    
    x_err(i) = norm(x_dif);
    v_err(i) = norm(v_dif);

    % Proportional, Integral, Derivative terms
    P = K_p * error;
    I = K_i * sum_error * dt;
    D = K_d * (error - prev_error) / dt;

    % Control input
    u = P + I + D;
%     u
    if norm(u) > u_limit
        disp(['step ' num2str(i) ' over limit'])
        u = u * u_limit / norm(u);
    end
%     u

    % Update sum of errors and previous error for the next iteration
    sum_error = sum_error + error;
    prev_error = error;
    
    if norm(x_dif) < dist_goal & norm(v_dif) < vel_goal
        T_true = i;
        break;
    end

end

%%

figure(fig_ind)
subplot(2, 3, 1)
hold on
scatter3(x1(:, 1), x1(:, 2), x1(:, 3), '.')
scatter3(x2(:, 1), x2(:, 2), x2(:, 3), '.')
title('dog chase')

subplot(2, 3, 2)
hold on
plot((1:T)/20, x_err)
xlim([0, T_true/fps])
ylim([0, 150])
xlabel('time [sec]')
ylabel('distance [m]')
title('Real Distance [m]')

subplot(2, 3, 3)
hold on
plot((1:T)/20, v_err)
xlim([0, T_true/fps])
%ylim([0, 200])
xlabel('time [sec]')
ylabel('velocity [m/sec]')
title('Real Velocity Difference')

subplot(2, 3, 4)
hold on
plot((1:T)/20, SNR_arr)
xlim([0, T_true/fps])
xlabel('time [sec]')
ylabel('SNR [dB]')
title('SNR')

subplot(2, 3, 5)
hold on
plot((1:T)/20, dist_err)
xlim([0, T_true/fps])
xlabel('time [sec]')
ylabel('distance error [m]')
title('Distance Estimation Error')

subplot(2, 3, 6)
hold on
plot((1:T)/20, DOA_err)
xlim([0, T_true/fps])
xlabel('time [sec]')
ylabel('3D angle error [theta]')
title('Angle Estimation Error')

%%

% if estim_dist_mode
%     if estim_DOA_mode
%         estim_text = '_estim dist & DOA';
%     else
%         estim_text = '_estim dist';
%     end
% else
%     if estim_DOA_mode
%         estim_text = '_estim DOA';
%     else
%         estim_text = 'true_loc';
%     end
% end
% 
% graph_file_name = ['tracking_' num2str(K_p) '_' num2str(K_i) '_' num2str(K_d) '_' estim_text '.mat'];
% save(fullfile("data/graphs/tracking", graph_file_name), "x1", "x2", "x_err", "v_err", "SNR_arr", "dist_err", "DOA_err");

%%
% for i = 1:T
%     
%     if i > 1
%         v1(i, :) = v1(i-1, :) + u;
%         x1(i, :) = x1(i-1, :) + v1(i, :);
%         x2(i, :) = x2(i-1, :) + v2(i, :);
%     end
%     
%     
%     ct = 1;
%     cr = 0.05;
%     
%     x_dif = x2(i, :) - x1(i, :);
%     v_dif = v1(i, :) - v2(i, :);
%     
%     alpha = acos(x_dif * v_dif.' / (norm(x_dif) * norm(v_dif)));
% %     beta = acos(x_dif(1) / norm(x_dif));
%     
%     u_rad = cr * tanh(log(norm(x_dif))) * (x_dif / norm(x_dif));
%     u_tang = ct * (-sin(alpha)^2) * [u_rad(2), -u_rad(1)];
%     u = u_tang + u_rad;
%     
% end

