function [estim_DOA, estim_DOA_ind] = calc_RobustMax(SRP, DOA, max_per, N)

    % [SRP_sorted, SRP_ind] = sort(SRP, 'descend');
    % DOA_sorted = DOA(SRP_ind, :);
    SRP_ind = find(SRP > max(SRP) * max_per);
    SRP_trimmed = SRP(SRP_ind, :);
    DOA_trimmed = DOA(SRP_ind, :);

%     figure(1);
%     scatter3(DOA(:, 1), DOA(:, 2), DOA(:, 3), 70, transpose(SRP), 'filled');axis equal;
%     cb = colorbar;
% 
%     figure(2);
%     scatter3(DOA_trimmed(:, 1), DOA_trimmed(:, 2), DOA_trimmed(:, 3), 70, transpose(SRP_trimmed), 'filled');axis equal;
%     cb = colorbar;

    % DOA_trimmed = DOA_sorted(1:RM_param, :);
    % N = 7;
    near_ind = zeros(N, size(SRP, 1));
    SRP_avg = zeros(size(SRP, 1), 1);

    for i = SRP_ind.'

        DOA_cent = DOA - repmat(DOA(i, :), size(DOA, 1), 1);
        [~, near_ind(:, i)] = mink(vecnorm(DOA_cent, 2, 2), N, 'ComparisonMethod', 'abs');
        SRP_avg(i) = mean(SRP(near_ind(:, i)));
    end

    DOA_trimmed_avg = DOA(find(SRP_avg), :);
%     figure(3);
%     scatter3(DOA_trimmed(:, 1), DOA_trimmed(:, 2), DOA_trimmed(:, 3), 70, transpose(SRP_avg(find(SRP_avg))), 'filled');axis equal;
%     cb = colorbar;

    [~, estim_DOA_ind] = max(SRP_avg);
    estim_DOA = DOA(estim_DOA_ind, :);

end

% [~, DOA_cent_ind] = min(norm(DOA_trimmed - mean(DOA_trimmed)));
% DOA_cent = DOA_trimmed(DOA_cent_ind, :);
% DOA_centralized = DOA_trimmed - DOA_cent;
% 
% DOA_clust_ind = find(vecnorm(DOA_centralized, 2, 2) < 0.5);
% 
% out_ind = zeros(size(DOA_clust_ind, 1), 1);
% 
% [~, max_ind] = max(SRP_trimmed(DOA_clust_ind));
% estim_DOA = DOA_trimmed(max_ind, :);

% for i = 1:size(DOA_clust_ind, 1)
%     X = (DOA(:, 1) == DOA_trimmed(DOA_clust_ind(i), 1));
%     Y = (DOA(:, 2) == DOA_trimmed(DOA_clust_ind(i), 2));
%     Z = (DOA(:, 3) == DOA_trimmed(DOA_clust_ind(i), 3));
%     out_ind(i) = find(X & Y & Z);
% end




% DOA_norm_sum = zeros(RM_param, 1);
% 
% for i = 1:RM_param
%     for j = 1:RM_param
%         DOA_norm_sum(i) = DOA_norm_sum(i) + norm(DOA_trimmed(i) - DOA_trimmed(j), 2);
%     end
% end
% 
% [~, min_ind] = min(DOA_norm_sum);
% [maxSRP, max_ind] = max(SRP);
% 
% if SRP_sorted(min_ind) * bar > maxSRP
%     estim_DOA = DOA_trimmed(min_ind, :);
% else
%     estim_DOA = DOA_trimmed(1, :);
% 
% end











% [SRP_sorted, SRP_ind] = sort(SRP, 'descend');
% DOA_sorted = DOA(SRP_ind, :);
% DOA_trimmed = DOA_sorted(1:RM_param, :);
% 
% DOA_norm_sum = zeros(RM_param, 1);
% 
% for i = 1:RM_param
%     for j = 1:RM_param
%         DOA_norm_sum(i) = DOA_norm_sum(i) + norm(DOA_trimmed(i) - DOA_trimmed(j), 2);
%     end
% end
% 
% [~, min_ind] = min(DOA_norm_sum);
% [maxSRP, max_ind] = max(SRP);
% 
% if SRP_sorted(min_ind) * bar > maxSRP
%     estim_DOA = DOA_trimmed(min_ind, :);
% else
%     estim_DOA = DOA_trimmed(1, :);
% 
% end




