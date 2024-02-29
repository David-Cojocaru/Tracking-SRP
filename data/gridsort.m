

[V, F] = icosphere(16);
[v1, f1] = icosphere(1);


J = zeros(length(V(:,1)),3,length(v1(:,1)));
N = floor(length(V(:,1))/12);
L = zeros(N,3,length(v1(:,1)));

for i = 1:length(v1)
    
    P = repmat(v1(i,:), size(V,1), 1);
    V_cent = V - P;
    
%     V_cent_sorted = sortrows(V_cent, 'ComparisonMethod','abs');
    [~, V_index] = sort(vecnorm(V_cent, 2, 2), 'ascend');
    J(:,:,i) = V_cent(V_index, :) + P;
    L(:,:,i) = J(1:N,:,i);
end

hold on
scatter3(L(:,1,1),L(:,2,1),L(:,3,1));axis equal;
trimesh(icosphere(16));axis equal; 