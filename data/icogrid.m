
% FINAL = 8
% B = {};
% for i = 0:FINAL
%     tic;
%     [B{end + 1}, ~] = icosphere(2^i);
%     toc;
% end
% 
% save(fullfile("data/icogrids", 'base_icogrid.mat'), "B");

%%

% figure();
% for i = 0:5
%     subplot(2, 3, i + 1)
%     title(['Degree' i])
%     trimesh(icosphere(2^i));axis equal;
% end

%%

% trimesh(icosphere(2^4));axis equal;% grid off
% base_angle_res = 63.4290;
% res_arr = [2^0, 2^1, 2^2, 2^3, 2^4, 2^5, 2^6, 2^7];
% angle_res = base_angle_res./res_arr;

 %%
 
tmp = load(fullfile("data/icogrids", ['base_icogrid.mat']));
B = tmp.B;

alpha = 0;

only_front_half = 1;
front_axis = 1;

for res_final = 6
    
    for res_start = 1%0:res_final
        tic;
        for res_num = 1%1:res_final - res_start + 1
            
            if res_num ~= 1
                if mod((res_final - res_start), (res_num - 1)) > 0
                    continue;
                end
            end
            
            if res_num == 1
                res_stages = res_final;
                res_jump = 1;
            else
                res_jump = max(1, floor((res_final - res_start)/(res_num - 1)));
                res_stages = res_start:res_jump:res_final;
            end

            % res_stages = 4;

            res_index = 2 .^ res_stages;
            tran_res = 2;
            M = {};
            G = {};

            index = 0;
%             figure();
            for n = res_stages

                index = index + 1;

                G{end + 1} = B{1, res_stages(index) + 1};

                if index == 1
                    
                    if alpha < 90
                        
                        L = B{1, res_stages(index) + 1};
                        L_sph = zeros(size(L));
                        [L_sph(:, 1), L_sph(:, 2), L_sph(:, 3)] = cart2sph(L(:, 1), L(:, 2), L(:, 3));
                        if only_front_half
                            L_ind = find(L_sph(:, 2) >= deg2rad(-alpha) & sign(front_axis)*L(:, abs(front_axis)) >= 0);
                        else
                            L_ind = find(L_sph(:, 2) >= deg2rad(-alpha));
                        end
                        
                        M{end + 1} = L(L_ind, :);
                        
%                         hold on;
%                         scatter3(L(L_ind,1,1),L(L_ind,2,1),L(L_ind,3,1));axis equal;
%                         trimesh(icosphere(res_index(index)));axis equal; 
 
                        
                    else
                        
                        M{end + 1} = B{1, res_stages(index) + 1};
                        
                    end

                else

                    V = B{1, res_stages(index) + 1};
                    v = B{1, res_stages(index - 1) + 1};
                    L = M{1, end};
%                     figure(index);
                    M{end + 1} = gridsort(V, v, L, res_index(index), res_index(index - 1), tran_res, res_jump);


                end

            end
            
            if res_num == 1
                if alpha < 90
                    if only_front_half
                        save(fullfile("data/icogrids", ['icogrid_' num2str(res_final) '_alpha_' num2str(alpha) '_front axis_' num2str(front_axis) '.mat']), "M");
                    else
                        save(fullfile("data/icogrids", ['icogrid_' num2str(res_final) '_alpha_' num2str(alpha) '.mat']), "M");
                    end
                else
                    if only_front_half
                        save(fullfile("data/icogrids", ['icogrid_' num2str(res_final) '_front axis_' num2str(front_axis) '.mat']), "M");
                    else
                        save(fullfile("data/icogrids", ['icogrid_' num2str(res_final) '.mat']), "M");
                    end
                end
                
%                 % grid_dim = icogrid
%                 tmp = load(fullfile("data/icogrids", ['icogrid_' num2str(res_final) '.mat']));
%                 % G = tmp.G;
%                 M = tmp.M;
% 
%                 tmp = load(fullfile("data/mic_arrays", "anti_prism_mic_array.mat"));
%                 micPos = tmp.micPos;
% 
%                 c = 343;
%                 DOA_list = {};
%                 Delta_list = {};
%                 m = size(micPos,1);
%                 P = m*(m-1)/2;
% 
%                 for i = 1:(size(M, 2))
% 
%                         L = M{1, i};
%                         DOA_i = L;
%                         Delta_t_i = zeros(size(L, 1), P, size(L, 3));
% 
%                         for j = 1:size(L, 3)
% 
%                                 Delta_t_i(:, :, j) = gen_searchIcoGrid(micPos, L(:, :, j), c);
%                         end
%                 %     end
% 
%                     Delta_list{end + 1} = Delta_t_i;
%                     DOA_list{end + 1} = DOA_i;
%                 end
% 
%                 save(fullfile("data/DOAs", ['DOAicogrid_' num2str(res_final) '.mat']), "DOA_list", "Delta_list");
                
            else
                if alpha < 90
                    if only_front_half
                        save(fullfile("data/icogrids", ['icogrid_' num2str(res_start) '_' num2str(res_num) '_' num2str(res_final) '_alpha_' num2str(alpha) '_front axis_' num2str(front_axis) '.mat']), "M");
                    else
                        save(fullfile("data/icogrids", ['icogrid_' num2str(res_start) '_' num2str(res_num) '_' num2str(res_final) '_alpha_' num2str(alpha) '.mat']), "M");
                    end
                else
                    if only_front_half
                        save(fullfile("data/icogrids", ['icogrid_' num2str(res_start) '_' num2str(res_num) '_' num2str(res_final) '_front axis_' num2str(front_axis) '.mat']), "M");
                    else
                        save(fullfile("data/icogrids", ['icogrid_' num2str(res_start) '_' num2str(res_num) '_' num2str(res_final) '.mat']), "M");
                    end
                end
                
%                 % grid_dim = icogrid
%                 tmp = load(fullfile("data/icogrids", ['icogrid_' num2str(res_start) '_' num2str(res_num) '_' num2str(res_final) '.mat']));
%                 % G = tmp.G;
%                 M = tmp.M;
% 
%                 tmp = load(fullfile("data/mic_arrays", "platform_mic_array.mat"));
%                 micPos = tmp.micPos;
% 
%                 c = 343;
%                 DOA_list = {};
%                 Delta_list = {};
%                 m = size(micPos,1);
%                 P = m*(m-1)/2;
% 
%                 for i = 1:(size(M, 2))
% 
%                         L = M{1, i};
%                         DOA_i = L;
%                         Delta_t_i = zeros(size(L, 1), P, size(L, 3));
% 
%                         for j = 1:size(L, 3)
% 
%                                 Delta_t_i(:, :, j) = gen_searchIcoGrid(micPos, L(:, :, j), c);
%                         end
%                 %     end
% 
%                     Delta_list{end + 1} = Delta_t_i;
%                     DOA_list{end + 1} = DOA_i;
%                 end
% 
%                 save(fullfile("data/DOAs", ['DOAicogrid_' num2str(res_start) '_' num2str(res_num) '_' num2str(res_final) '.mat']), "DOA_list", "Delta_list");
            end
             
        end
        toc;
    end
    
end

%%



%%

function L = gridsort(V, v, l, m, n, tran_res, res_jump)

    
%     V = B{1, k};
%     v = B{1, l};  
%     size(V, 1)
%     size(l, 1) * size(l, 3)
    %J = zeros(size(V, 1), 3, size(l, 1) * size(l, 3));
    N = 3*(res_jump + 1)^2 - 3*(res_jump + 1) + 1; %floor(size(V, 1) * tran_res / size(v, 1));
    L = zeros(N, 3, size(l, 1) * size(l, 3));
    
    k = 0;

    for i = 1:size(l, 3)
        
        for j = 1:size(l, 1)
            
            k = k + 1;

            P = repmat(l(j,:, i), size(V,1), 1);
            V_cent = V - P;

            [~, V_index] = mink(vecnorm(V_cent, 2, 2), N);
            L(:,:,k) = V_cent(V_index, :) + P(V_index, :);
          
        end
    end

%     for i = 1:size(v, 1)
% 
%         P = repmat(v(i,:), size(V,1), 1);
%         V_cent = V - P;
% 
%         [~, V_index] = sort(vecnorm(V_cent, 2, 2), 'ascend');
%         J(:,:,i) = V_cent(V_index, :) + P;
%         L(:,:,i) = J(1:N,:,i);
%           
%     end

%     k = 0;
% 
%     for i = 1:size(l, 3)
%         
%         for j = 1:size(l, 1)
%             
%             k = k + 1;
% 
%             P = repmat(l(j,:, i), size(V,1), 1);
%             V_cent = V - P;
% 
%             [~, V_index] = sort(vecnorm(V_cent, 2, 2), 'ascend');
%             J(:,:,k) = V_cent(V_index, :) + P;
%             L(:,:,k) = J(1:N,:,k);
%           
%         end
%     end

    
%     hold on
%     for i = 1:size(L, 3)
%         scatter3(L(:,1,i),L(:,2,i),L(:,3,i));axis equal;
%     end
%     scatter3(L(:,1,1),L(:,2,1),L(:,3,1));axis equal;
%     trimesh(icosphere(m));axis equal; 
 
end

    
    