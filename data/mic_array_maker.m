
c = 343;
%% anti-prism ensamble:
anti_prism_height = 0.1;
micPos = [];
micPos = [micPos; [0.1, 0.1, 0]];
micPos = [micPos; [-0.1, 0.1, 0]];
micPos = [micPos; [-0.1, -0.1, 0]];
micPos = [micPos; [0.1, -0.1, 0]];
micPos = [micPos; [0.05, 0, anti_prism_height]];
micPos = [micPos; [0, 0.05, anti_prism_height]];
micPos = [micPos; [-0.05, 0, anti_prism_height]];
micPos = [micPos; [0, -0.05, anti_prism_height]];
micPos = micPos - [0, 0, anti_prism_height/2];
% micPos = micPos*45;

arrayCenterPos = sum(micPos, 1) / size(micPos, 1);

dist_mat = zeros(size(micPos, 1), size(micPos, 1));
for m = 1:size(micPos, 1)
    for mprime = 1:size(micPos, 1)
        dist_mat(m, mprime) = norm(micPos(m, :) - micPos(mprime, :));
    end
end

save(fullfile("data/mic_arrays", "anti_prism_mic_array.mat"), "micPos", "arrayCenterPos", "dist_mat");

%% anti-prism-double-elevation ensamble:
micPos = [];
micPos = [micPos; [0.1, 0.1, 0]];
micPos = [micPos; [-0.1, 0.1, 0]];
micPos = [micPos; [-0.1, -0.1, 0]];
micPos = [micPos; [0.1, -0.1, 0]];
micPos = [micPos; [0.05, 0, anti_prism_height * 2]];
micPos = [micPos; [0, 0.05, anti_prism_height]];
micPos = [micPos; [-0.05, 0, anti_prism_height * 2]];
micPos = [micPos; [0, -0.05, anti_prism_height]];
micPos = micPos - [0, 0, 0.05];

arrayCenterPos = sum(micPos, 1) / size(micPos, 1);

dist_mat = zeros(size(micPos, 1), size(micPos, 1));
for m = 1:size(micPos, 1)
    for mprime = 1:size(micPos, 1)
        dist_mat(m, mprime) = norm(micPos(m, :) - micPos(mprime, :));
    end
end

save(fullfile("data/mic_arrays", "anti_prism_2_mic_array.mat"), "micPos", "arrayCenterPos", "dist_mat");

%% cube ensamble

cube_side = 1;
micPos = [];
micPos = [micPos; [0, 0, 0]];
micPos = [micPos; [cube_side, 0, 0]];
micPos = [micPos; [0, cube_side, 0]];
micPos = [micPos; [0, 0, cube_side]];
micPos = [micPos; [cube_side, cube_side, 0]];
micPos = [micPos; [cube_side, 0, cube_side]];
micPos = [micPos; [0, cube_side, cube_side]];
micPos = [micPos; [cube_side, cube_side, cube_side]];
micPos = micPos - sum(micPos, 1) / size(micPos, 1);

arrayCenterPos = sum(micPos, 1) / size(micPos, 1);

dist_mat = zeros(size(micPos, 1), size(micPos, 1));
for m = 1:size(micPos, 1)
    for mprime = 1:size(micPos, 1)
        dist_mat(m, mprime) = norm(micPos(m, :) - micPos(mprime, :));
    end
end

save(fullfile("data/mic_arrays", "cube_mic_array.mat"), "micPos", "arrayCenterPos", "dist_mat");

%% hexadome ensamble

hexadome_height = 0.02;
micPos = [];
micPos = [micPos; [0.2, 0, 0]];
micPos = [micPos; [0.1, 0.1732, 0]];
micPos = [micPos; [-0.1, 0.1732, 0]];
micPos = [micPos; [-0.2, 0, 0]];
micPos = [micPos; [-0.1, -0.1732, 0]];
micPos = [micPos; [0.1, -0.1732, 0]];
micPos = [micPos; [0.1, 0, hexadome_height]];
micPos = [micPos; [-0.1, 0, hexadome_height]];
micPos(:, 3) = micPos(:, 3) - hexadome_height/4; 

arrayCenterPos = sum(micPos, 1) / size(micPos, 1);
% [DOAvec_i, Delta_t_i] = gen_searchGrid(micPos, ang_pol, ang_az, 'spherical', c);

dist_mat = zeros(size(micPos, 1), size(micPos, 1));
for m = 1:size(micPos, 1)
    for mprime = 1:size(micPos, 1)
        dist_mat(m, mprime) = norm(micPos(m, :) - micPos(mprime, :));
    end
end

save(fullfile("data/mic_arrays", "hexadome_mic_array.mat"), "micPos", "arrayCenterPos", "dist_mat");

%% billiardome ensamble

billiardome_height = 0.02;
micPos = [];
micPos = [micPos; [-0.1, 0.05, 0]];
micPos = [micPos; [0, 0.05, 0]];
micPos = [micPos; [0.1, 0.05, 0]];
micPos = [micPos; [-0.1, -0.05, 0]];
micPos = [micPos; [0, -0.05, 0]];
micPos = [micPos; [0.1, -0.05, 0]];
micPos = [micPos; [0.05, 0, billiardome_height]];
micPos = [micPos; [-0.05, 0, billiardome_height]];
micPos(:, 3) = micPos(:, 3) - billiardome_height/4; 

arrayCenterPos = sum(micPos, 1) / size(micPos, 1);
% [DOAvec_i, Delta_t_i] = gen_searchGrid(micPos, ang_pol, ang_az, 'spherical', c);

dist_mat = zeros(size(micPos, 1), size(micPos, 1));
for m = 1:size(micPos, 1)
    for mprime = 1:size(micPos, 1)
        dist_mat(m, mprime) = norm(micPos(m, :) - micPos(mprime, :));
    end
end

save(fullfile("data/mic_arrays", "billiardome_mic_array.mat"), "micPos", "arrayCenterPos", "dist_mat");

%% platform ensamble

platform_height = 0.02;
micPos = [];
micPos = [micPos; [0.3, 0.05, 0]];
micPos = [micPos; [0, 0.055, 0]];
micPos = [micPos; [-0.25, 0.055, 0]];
micPos = [micPos; [0.3, -0.05, 0]];
micPos = [micPos; [0, -0.055, 0]];
micPos = [micPos; [-0.25, -0.055, 0]];
micPos = [micPos; [0, 0, platform_height]];
micPos = [micPos; [0.3, 0, platform_height]]; 

arrayCenterPos = sum(micPos, 1) / size(micPos, 1);
% [DOAvec_i, Delta_t_i] = gen_searchGrid(micPos, ang_pol, ang_az, 'spherical', c);

dist_mat = zeros(size(micPos, 1), size(micPos, 1));
for m = 1:size(micPos, 1)
    for mprime = 1:size(micPos, 1)
        dist_mat(m, mprime) = norm(micPos(m, :) - micPos(mprime, :));
    end
end

save(fullfile("data/mic_arrays", "platform_mic_array.mat"), "micPos", "arrayCenterPos", "dist_mat");

%% tube ensamble

tube_rad = 0.1;
tube_len = 0.7;
mic_offset = (tube_len - 0.1) / 2;
tube_height = tube_rad - 0.03;
micPos = [];
micPos = [micPos; [mic_offset, tube_rad, 0]];
micPos = [micPos; [0, tube_rad, 0]];
micPos = [micPos; [-mic_offset, tube_rad, 0]];
micPos = [micPos; [mic_offset, -tube_rad, 0]];
micPos = [micPos; [0, -tube_rad, 0]];
micPos = [micPos; [-mic_offset, -tube_rad, 0]];
micPos = [micPos; [tube_len/2, 0, tube_height]];
micPos = [micPos; [-tube_len/2, 0, tube_height]]; 

arrayCenterPos = sum(micPos, 1) / size(micPos, 1);
% [DOAvec_i, Delta_t_i] = gen_searchGrid(micPos, ang_pol, ang_az, 'spherical', c);

dist_mat = zeros(size(micPos, 1), size(micPos, 1));
for m = 1:size(micPos, 1)
    for mprime = 1:size(micPos, 1)
        dist_mat(m, mprime) = norm(micPos(m, :) - micPos(mprime, :));
    end
end

save(fullfile("data/mic_arrays", "tube_mic_array.mat"), "micPos", "arrayCenterPos", "dist_mat");

%% octahedron ensamble

octahedron_height = 0.1;
micPos = [];
micPos = [micPos; [-0.1, 0, 0]];
micPos = [micPos; [0, 0.1, 0]];
micPos = [micPos; [0.1, 0, 0]];
micPos = [micPos; [0, -0.1, 0]];
micPos = [micPos; [0, 0, octahedron_height]];
micPos = [micPos; [0, 0, -octahedron_height]];

arrayCenterPos = sum(micPos, 1) / size(micPos, 1);
% [DOAvec_i, Delta_t_i] = gen_searchGrid(micPos, ang_pol, ang_az, 'spherical', c);

dist_mat = zeros(size(micPos, 1), size(micPos, 1));
for m = 1:size(micPos, 1)
    for mprime = 1:size(micPos, 1)
        dist_mat(m, mprime) = norm(micPos(m, :) - micPos(mprime, :));
    end
end

save(fullfile("data/mic_arrays", "octahedron_mic_array.mat"), "micPos", "arrayCenterPos", "dist_mat");

%% tetrahedron ensamble

tetrahedron_side = 0.1;
micPos = [];
micPos = [micPos; [0, 0, 0]];
micPos = [micPos; [tetrahedron_side, tetrahedron_side, 0]];
micPos = [micPos; [tetrahedron_side, 0, tetrahedron_side]];
micPos = [micPos; [0, tetrahedron_side, tetrahedron_side]];
micPos = micPos / sqrt(2);

micPos = micPos - sum(micPos, 1) / size(micPos, 1);
% scatter3(micPos(:,1),micPos(:,2),micPos(:,3));axis equal;

arrayCenterPos = sum(micPos, 1) / size(micPos, 1);
% [DOAvec_i, Delta_t_i] = gen_searchGrid(micPos, ang_pol, ang_az, 'spherical', c);

save(fullfile("data/mic_arrays", "tetrahedron_mic_array.mat"), "micPos", "arrayCenterPos", "dist_mat");
%% icosahedron ensamble

[micPos, ~] = icosphere(1);

micPos = micPos./10;

arrayCenterPos = sum(micPos, 1) / size(micPos, 1);
% [DOAvec_i, Delta_t_i] = gen_searchGrid(micPos, ang_pol, ang_az, 'spherical', c);

dist_mat = zeros(size(micPos, 1), size(micPos, 1));
for m = 1:size(micPos, 1)
    for mprime = 1:size(micPos, 1)
        dist_mat(m, mprime) = norm(micPos(m, :) - micPos(mprime, :));
    end
end

save(fullfile("data/mic_arrays", "icosahedron.mat"), "micPos", "arrayCenterPos", "dist_mat");

%% Line

line_len = 1;
mic_offset = (line_len) / 7;

micPos = [];
micPos = [micPos; [0, 0, 0]];
micPos = [micPos; [mic_offset * 1, 0, 0]];
micPos = [micPos; [mic_offset * 2, 0, 0]];
micPos = [micPos; [mic_offset * 3, 0, 0]];
micPos = [micPos; [mic_offset * 4, 0, 0]];
micPos = [micPos; [mic_offset * 5, 0, 0]];
micPos = [micPos; [mic_offset * 6, 0, 0]];
micPos = [micPos; [mic_offset * 7, 0, 0]];
micPos(:, 1) = micPos(:, 1) - line_len/2; 

arrayCenterPos = sum(micPos, 1) / size(micPos, 1);
% [DOAvec_i, Delta_t_i] = gen_searchGrid(micPos, ang_pol, ang_az, 'spherical', c);

dist_mat = zeros(size(micPos, 1), size(micPos, 1));
for m = 1:size(micPos, 1)
    for mprime = 1:size(micPos, 1)
        dist_mat(m, mprime) = norm(micPos(m, :) - micPos(mprime, :));
    end
end

save(fullfile("data/mic_arrays", "line_mic_array.mat"), "micPos", "arrayCenterPos", "dist_mat");


%% Circle 

% circ_rad = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3, 3.5, 4];
circ_rad = [0.3, 0.5, 0.75];

for i = 1:size(circ_rad, 2)
    theta = 2*pi / 4;
    theta_off = pi / 20;
    r = ones(8, 1) * circ_rad(i);
    theta_mics = [theta * 0 - theta_off, theta * 0 + theta_off, theta * 1 - theta_off, theta * 1 + theta_off, theta * 2 - theta_off, theta * 2 + theta_off, theta * 3 - theta_off, theta * 3 + theta_off]'
    z = zeros(8, 1);

    % micPos = [];
    % micPos = [micPos; [1, 0, 0]];
    % micPos = [micPos; [pol2cart(theta * 1, 1), 0]];
    % micPos = [micPos; [pol2cart(theta * 2, 1), 0]];
    % micPos = [micPos; [pol2cart(theta * 3, 1), 0]];
    % micPos = [micPos; [pol2cart(theta * 4, 1), 0]];
    % micPos = [micPos; [pol2cart(theta * 5, 1), 0]];
    % micPos = [micPos; [pol2cart(theta * 6, 1), 0]];
    % micPos = [micPos; [pol2cart(theta * 7, 1), 0]];
    [x, y, z] = pol2cart(theta_mics, r, z);
    micPos = [x, y, z];
%     scatter3(micPos(:,1),micPos(:,2),micPos(:,3));axis equal;

    arrayCenterPos = sum(micPos, 1) / size(micPos, 1);
    % [DOAvec_i, Delta_t_i] = gen_searchGrid(micPos, ang_pol, ang_az, 'spherical', c);
    
    dist_mat = zeros(size(micPos, 1), size(micPos, 1));
    for m = 1:size(micPos, 1)
        for mprime = 1:size(micPos, 1)
            dist_mat(m, mprime) = norm(micPos(m, :) - micPos(mprime, :));
        end
    end

    save(fullfile("data/mic_arrays", ['circle_' num2str(circ_rad(i)) '_pairs_mic_array.mat']), "micPos", "arrayCenterPos", "dist_mat");
end

%% H-323

base_len = 1;
height = 0.5;
mic_offset_b = base_len / 2;
mic_offset_h = height / 3;

micPos = [];
micPos = [micPos; [mic_offset_b, mic_offset_h * 3, 0]];
micPos = [micPos; [0, mic_offset_h * 3, 0]];
micPos = [micPos; [-mic_offset_b, mic_offset_h * 3, 0]];
micPos = [micPos; [0, mic_offset_h * 2, 0]];
micPos = [micPos; [0, mic_offset_h * 1, 0]];
micPos = [micPos; [mic_offset_b, 0, 0]];
micPos = [micPos; [0, 0, 0]];
micPos = [micPos; [-mic_offset_b, 0, 0]];
micPos(:, 2) = micPos(:, 2) - height/2;
% micPos = micPos/4;
% scatter3(micPos(:,1),micPos(:,2),micPos(:,3));axis equal;

arrayCenterPos = sum(micPos, 1) / size(micPos, 1);
% [DOAvec_i, Delta_t_i] = gen_searchGrid(micPos, ang_pol, ang_az, 'spherical', c);

dist_mat = zeros(size(micPos, 1), size(micPos, 1));
for m = 1:size(micPos, 1)
    for mprime = 1:size(micPos, 1)
        dist_mat(m, mprime) = norm(micPos(m, :) - micPos(mprime, :));
    end
end

save(fullfile("data/mic_arrays", "H323_mic_array.mat"), "micPos", "arrayCenterPos", "dist_mat");

%% H-242

base_len = 2;
height = 1;
mic_offset_b = base_len / 2;
mic_offset_h = height / 5;

micPos = [];
micPos = [micPos; [mic_offset_b, mic_offset_h * 5, 0]];
micPos = [micPos; [0, mic_offset_h * 4, 0]];
micPos = [micPos; [-mic_offset_b, mic_offset_h * 5, 0]];
micPos = [micPos; [0, mic_offset_h * 3, 0]];
micPos = [micPos; [0, mic_offset_h * 2, 0]];
micPos = [micPos; [mic_offset_b, 0, 0]];
micPos = [micPos; [0, mic_offset_h * 1, 0]];
micPos = [micPos; [-mic_offset_b, 0, 0]];
micPos(:, 2) = micPos(:, 2) - height/2;
% scatter3(micPos(:,1),micPos(:,2),micPos(:,3));axis equal;

arrayCenterPos = sum(micPos, 1) / size(micPos, 1);
% [DOAvec_i, Delta_t_i] = gen_searchGrid(micPos, ang_pol, ang_az, 'spherical', c);

dist_mat = zeros(size(micPos, 1), size(micPos, 1));
for m = 1:size(micPos, 1)
    for mprime = 1:size(micPos, 1)
        dist_mat(m, mprime) = norm(micPos(m, :) - micPos(mprime, :));
    end
end

save(fullfile("data/mic_arrays", "H242_mic_array.mat"), "micPos", "arrayCenterPos", "dist_mat");


