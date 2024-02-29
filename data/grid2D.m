% c = 343;
% 
% tmp = load(fullfile("data/mic_arrays", "platform_mic_array.mat"));
% micPos = tmp.micPos;
% 
% B = {};
% 
% for grid_res = 1:90
%     
%     if mod(360, grid_res) == 0
%         % polar angles of candidate locations
%         ang_pol = 90;
%         % azimuth angles of candidate locations 
%         ang_az = 0:grid_res:359;
%         % compute candidate DOA vectors and TDOAs
%         % [DOAvec_i, Delta_t_i] = gen_searchGrid(micPos, ang_pol, ang_az, 'spherical', c);
% 
%         DOA_vec = zeros(size(ang_az, 2), 3);
% 
%         i = 0;
% 
%         for n_dim1 = ang_pol
%             for n_dim2 = ang_az
% 
%                 i = i + 1;
% 
%                 DOA_vec(i, :) = [sin(deg2rad(n_dim1))*cos(deg2rad(n_dim2)),...
%                 sin(deg2rad(n_dim1))*sin(deg2rad(n_dim2)),...
%                 cos(deg2rad(n_dim1))];
% 
%             end
%         end
%         
%         B{1, end+1} = DOA_vec;
%     end
% end
% 
% save(fullfile("data/grids2D", 'base_grids.mat'), "B");

%%

tmp = load(fullfile("data/grids2D", 'base_grids.mat'));
B = tmp.B;
final = 1;
for first = 1:size(B, 2)
    M = {};
    res_stages = [first, final];
    
    if first == final
        res_num = 1;
    else
        res_num = size(res_stages, 2);
    end

    for n = 1:res_num

        res_fine = size(B{1, res_stages(n)}, 1);

        if n == 1

            M{end + 1} = B{1, res_stages(n)};

        else

            res_coarse = size(B{1, res_stages(n-1)}, 1);

            V = B{1, res_stages(n)};

            l = M{1, end};
            N = floor(res_fine/res_coarse);
            L = zeros(N, 3, size(l, 1) * size(l, 3));

            k = 0;

            for i = 1:size(l, 3)

                for j = 1:size(l, 1)

                    k = k + 1;

                    P = repmat(l(j,:, i), size(V, 1), 1);
                    V_cent = V - P;

                    [~, V_index] = mink(vecnorm(V_cent, 2, 2), N);
                    L(:,:,k) = V_cent(V_index, :) + P(V_index, :);

                end
            end

            M{end + 1} = L;

        end

    end

    if res_num == 1
        save(fullfile("data/grids2D", ['grid_' num2str(size(B{1, res_stages(end)}, 1)) '.mat']), "M");

        tmp = load(fullfile("data/grids2D", "grid_24_2_360.mat"));
        % G = tmp.G;
        M = tmp.M;

        tmp = load(fullfile("data/mic_arrays", "platform_mic_array.mat"));
        micPos = tmp.micPos;

        c = 343;
        DOA_list = {};
        Delta_list = {};
        m = size(micPos,1);
        P = m*(m-1)/2;

        for i = 1:(size(M, 2))

                L = M{1, i};
                DOA_i = L;
                Delta_t_i = zeros(size(L, 1), P, size(L, 3));

                for j = 1:size(L, 3)

                        Delta_t_i(:, :, j) = gen_searchIcoGrid(micPos, L(:, :, j), c);
                end

            Delta_list{end + 1} = Delta_t_i;
            DOA_list{end + 1} = DOA_i;
        end

        save(fullfile("data/DOAs", ['DOAgrid_' num2str(size(B{1, res_stages}, 1)) '.mat']), "DOA_list", "Delta_list");

    else
        save(fullfile("data/grids2D", ['grid_' num2str(size(B{1, res_stages(1)}, 1)) '_' num2str(res_num) '_' num2str(size(B{1, res_stages(end)}, 1)) '.mat']), "M");

        tmp = load(fullfile("data/grids2D", ['grid_' num2str(size(B{1, res_stages(1)}, 1)) '_' num2str(res_num) '_' num2str(size(B{1, res_stages(end)}, 1)) '.mat']));
        % G = tmp.G;
        M = tmp.M;

        tmp = load(fullfile("data/mic_arrays", "platform_mic_array.mat"));
        micPos = tmp.micPos;

        c = 343;
        DOA_list = {};
        Delta_list = {};
        m = size(micPos,1);
        P = m*(m-1)/2;

        for i = 1:(size(M, 2))

                L = M{1, i};
                DOA_i = L;
                Delta_t_i = zeros(size(L, 1), P, size(L, 3));

                for j = 1:size(L, 3)

                        Delta_t_i(:, :, j) = gen_searchIcoGrid(micPos, L(:, :, j), c);
                end

            Delta_list{end + 1} = Delta_t_i;
            DOA_list{end + 1} = DOA_i;
        end

        save(fullfile("data/DOAs", ['DOAgrid_' num2str(size(B{1, res_stages(1)}, 1)) '_' num2str(res_num) '_' num2str(size(B{1, res_stages(end)}, 1)) '.mat']), "DOA_list", "Delta_list");
    end
end


%%