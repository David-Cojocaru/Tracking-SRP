function [delta_t_i] = gen_searchIcoGrid(micPos, grid, c)

M = size(micPos,1);
P = M*(M-1)/2;
delta_t_i = zeros(size(grid, 1), P);

for i = 1:size(grid, 1)

    DOA_vec = grid(i, :);

    b = -transpose(DOA_vec);

    delta_t = zeros(P,1);
    p = 0;
    for mprime = 1:M
        for m = mprime+1:M
            p = p+1;

            a = micPos(m,:).' - micPos(mprime,:).';

            delta_t(p) = (norm(a)/c)*(a.'*b/(norm(a)*norm(b)));

        end
    end
    
    
    delta_t_i(i, :) = delta_t;
end

end