function [] = DOA_Delta_icogrid(mic_config, grid)

    tmp = load(fullfile("data/grids2D", [grid '.mat']));
    % G = tmp.G;
    M = tmp.M;

    tmp = load(fullfile("data/mic_arrays", "circle_mic_array.mat"));
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

    save(fullfile("data/DOAs", "DOA.mat"), "DOA_list", "Delta_list");

end
