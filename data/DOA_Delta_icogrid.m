function [] = DOA_Delta_icogrid(mic_config, icogrid)
    
    
    
    tmp = load(fullfile("data/icogrids", [icogrid '.mat']));
    M = tmp.M;

    tmp = load(fullfile("data/mic_arrays", [mic_config '.mat']));
    micPos = tmp.micPos;

    c = 343;
    m = size(micPos,1);
    P = m*(m-1)/2;
    hierarchy_depth = size(M, 2);
    Delta_list = cell(hierarchy_depth);
    DOA_list = cell(hierarchy_depth);

    for i = 1:(size(M, 2))

            L = M{1, i};
            DOA_i = L;
            Delta_t_i = zeros(size(L, 1), P, size(L, 3));

            for j = 1:size(L, 3)

                    Delta_t_i(:, :, j) = gen_searchIcoGrid(micPos, L(:, :, j), c);
            end

        Delta_list{i} = Delta_t_i;
        DOA_list{i} = DOA_i;
    end

    save(fullfile("data/DOAs", 'DOA.mat'), "DOA_list", "Delta_list");

end