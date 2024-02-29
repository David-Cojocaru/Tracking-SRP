%% PREAMBLE

clear;
cd(fileparts(mfilename('fullpath')));
addpath(genpath(pwd));

% run("icogrid.m");
% run("randPoints.m")
% run("mic_array_maker.m")
% run("DOA_Delta.m")


%% CONFIGURATION


%%% MODES/FLAGS

Attenuation_Mode = 1; % if enabled, the mono-channel input is passed through the attenuation filter to produce x_TD at distance r. The noise is callibrated to the amplitude of x_TD at distance r0, and can be adjusted via the Noise-Multiplier.
% if disabled, as long as the far-field assumption holds, the distance of the target does not affect anything, and the SNR can be set by hand.

DOA_processing = 1; % if enabled, the run will include processing of a new DOA file with a desired microphone configuration and search grid.

location_suffling = 1; % if enabled, the target locations will be randomly regenerated.

noised_input_mode = 1; % if enabled, the multi-channel input will be composed of the original signal and adative noise, with random delays between channels.

single_target = 1; % if enabled, the algorithm will run on a single target, recommended option for displaying intensity map for analysis along with enabling figure mode.

figure_mode = 1; % if enabled, every intensity map will be displayed on a figure.

search_grid_3D_mode = 1; % if enabled, the search grid will contain dots on the unit sphere, and if disabeld, on the xy unit circle. for the 2D mode we recommend disabeling the prediction 3D mode as well.

prediction_3D_mode = 1; % if enabled, the error will be messured as the angle between 2 3D-vectors, if disabled, the error will be messured as the angle between the 2D-vector projections (the azimuth component).

graph_mode = 0; % if enabled, for each target the error, SNR, and distance will be recorded and saved as a mat file, and can be plotted as a graph through scale_2_error.m

phat_mode = 0; % if enabled, the SRP algorithm uses PHAT transform as generelized weighting for the GCC

weight_gcc_mode = 0; % if enabled, the SRP algorithm uses distance-dependent heuristic as generelized weighting for the GCC


%%% COSTUM VARIABELS

mic_config = 'platform_mic_array'; % all microphone configuration files are found in the mic_array folder.

icogrid = 'icogrid_6_alpha_0_front axis_1';

grid_2D = 'grid_360';

signal_type = 'explosion';

noise_type = 'white_noise';

num_of_targets = 100;
furthest_distance_generated = 1000;
closest_distance_generated = 10;



%%% ACOUSTIC SETUP
% speed of sound
c = 343;
% sample rate
fs = 32000;
% bandlimit
percentage_of_pass_band = 32 / 32;
w_0 = pi*fs * percentage_of_pass_band;


if Attenuation_Mode
    noise_mul = -20; % Noise-Multiplier in dB
else
    SNR = 0; % SNR in dB
end


%%% MICROPHONE ARRAY
tmp = load(fullfile("data/mic_arrays", [mic_config '.mat']));
% array center
arrayCenterPos = tmp.arrayCenterPos;
dist_mat = tmp.dist_mat;
% microphone positions
micPos = tmp.micPos;
% number of microphones
M = size(micPos,1);
% number of microphones pairs
P = M*(M-1)/2;


%%% DOA LOAD

if DOA_processing
    if search_grid_3D_mode
        DOA_Delta_icogrid(mic_config, icogrid);
    else
        DOA_Delta_grid(mic_config, icogrid);
    end
end

tmp = load(fullfile("data/DOAs", 'DOA.mat')); % all DOA files and are found in the DOAs folder.
% the search grid itself
DOAvec_list = tmp.DOA_list;
% the matching pair-delays 
Delta_t_list = tmp.Delta_list;


%%% SOURCE LOCATIONS

if location_suffling
    randPoints(num_of_targets, furthest_distance_generated, closest_distance_generated);
end

tmp = load('randLoc.mat');
true_loc = tmp.true_loc;
% compute ground truth DOA vectors for source locations
true_DOAvec = calc_DOA(true_loc, arrayCenterPos);


%%% STFT PARAMETERS

win_time = 0.05; % 50ms
% window size
N_STFT = fs * win_time;
% shift
R_STFT = N_STFT/2;
% window
win = sqrt(hann(N_STFT,'periodic'));
N_STFT_half = floor(N_STFT/2)+1;
% frequency vector
omega = 2*pi*linspace(0,fs/2,N_STFT_half).';
% number of processed frames per location
L = 1;


%%% SRP APPROXIMATION PARAMETERS

% compute sampling period and number of samples within TDOA interval
[ T, N_mm ] = calc_sampleParam(micPos, w_0, c);
% number of auxilary samples
N_aux = 2;


%% PROCESSING

if single_target
    loc_num = 1;
else
    loc_num = size(true_loc,1);
end

% init results (per source location, frame, number of auxilary samples)
% approximation error in dB
res.approxErr_dB = zeros(loc_num, L, length(N_aux));
% localization error 
res.locErr = zeros(loc_num, L, length(N_aux)); % res.locErr(:,:,1) refers to  conventional SRP
res.locErrRM = zeros(loc_num, L, length(N_aux)); % res.locErr(:,:,1) refers to  conventional SRP


if graph_mode
    platform_err = zeros(loc_num, 1);
    platform_errRM = zeros(loc_num, 1);
    platform_SNR = zeros(loc_num, 1);
    platform_dist = zeros(loc_num, 1);
end


Tictoc = zeros(loc_num, 1);


for true_loc_idx = 1:loc_num

    
    disp(['PROCESSING SOURCE LOCATION ' num2str(true_loc_idx)])
    %%% GENERATE MICROPHONE SIGNALS    

    
    % Input Component
    
    [x_TDmono, Fs] = audioread([signal_type '.wav']); % The input signal as a stereo wav file

    x_TDmono = (x_TDmono(:, 1) + x_TDmono(:, 2)) / 2; % Converting the input to mono-channel wav file.

    fs_resample = fs; % The desired resampling frequency

    [x_TD, delay] = calc_INPUT_SIGNAL(x_TDmono, micPos, true_loc(true_loc_idx,:), c, Fs, fs_resample);
    
    multi_target_input_mode = 0;
    
    if multi_target_input_mode
        [x_TD_2, delay2] = calc_INPUT_SIGNAL(x_TDmono, micPos, true_loc(true_loc_idx + 1,:), c, Fs, fs_new);
        x_TD = x_TD + x_TD_2;
    end

    
    % noise component
    
    if noised_input_mode

        switch noise_type

            case 'white_noise'

                v_TD = randn(size(x_TD, 1) , size(x_TD, 2));

            case 'two_stroke'

                [v_TD, fsV] = audioread('two_stroke_noise.wav');
                v_TD = resample(v_TD, fs_resample, fsV);
                v_TD = repmat(v_TD(:, 1), 3, M);

                max_offset = 1000;
                noise_offset = floor(max_offset*rand(1,M));
                v_TD = delayseq(v_TD, noise_offset);
                v_TD = v_TD(1 + max_offset:size(x_TD, 1) + max_offset, :);

            case 'phantom'

                [v_TD, fsV] = audioread('phantom_4.wav');
                v_TD = resample(v_TD, fs_resample, fsV);
                v_TD = repmat(v_TD(:, 1), 3, M);

                max_offset = 1000;
                noise_offset = floor(max_offset*rand(1,M));
                v_TD = delayseq(v_TD, noise_offset);
                v_TD = v_TD(1 + max_offset:size(x_TD, 1) + max_offset, :);

        end
        
    end


    

    tic; % For each target the calculation time includes the STFT of the signal as well as the SRP algorithm, the simulation is not included since in real time usage of the algorithm, x_TD is a direct input.
    
    
    % transform to STFT domain
    
    x_STFT  = calc_STFT(x_TD, fs, win, N_STFT, R_STFT, 'onesided');
    
    if noised_input_mode
        
        v_STFT  = calc_STFT(v_TD, fs, win, N_STFT, R_STFT, 'onesided');
        
        if Attenuation_Mode
        
            dist_at_SNR_0 = 100;
            humidity = 50;

            [xr_STFT, SNR_prescale] = Attenuation(x_STFT, fs, norm(true_loc(true_loc_idx, :)), humidity, dist_at_SNR_0);
            [xr_STFT, v_STFT, scaling] = set_SNR(xr_STFT, v_STFT, SNR_prescale);

            v_STFT = db2mag(noise_mul) * v_STFT;
            SNR = SNR_prescale - noise_mul;
            disp(['SNR ' num2str(SNR)])
            
            if graph_mode
                graph_SNR(true_loc_idx) = SNR;
                graph_dist(true_loc_idx) = norm(true_loc(true_loc_idx, :));
            end

        else

            set_SNR(x_STFT, v_STFT, SNR);

        end
        
        
        % discard frames that do not contain real signal energy (local SNR below SNR barrier)
        SNR_barrier = SNR - 10;
        l = 1;
        useframe_idx = [];

        while length(useframe_idx) < L
            SNR_local = pow2db(sum(abs(xr_STFT(:,l,1)).^2)/sum(abs(v_STFT(:,l,1)).^2));
            if SNR_local > SNR_barrier
                useframe_idx = [useframe_idx, l];
            end
            l = l + 1;
        end
        
        
        % final microphone signal in STFT domain
        
        if Attenuation_Mode
            
            y_STFT = xr_STFT(:,useframe_idx,:) + v_STFT(:,useframe_idx,:);
            
        else
            
            y_STFT = x_STFT(:,useframe_idx,:) + v_STFT(:,useframe_idx,:);
            
        end
        
    else
        
        if Attenuation_Mode
            
            [xr_STFT, SNR_prescale] = Attenuation(x_STFT, fs, norm(true_loc(true_loc_idx, :)), humidity, dist_at_SNR_0);
            
            % final microphone signal in STFT domain
            y_STFT = xr_STFT(:,1,:)
            
        else
            
            % final microphone signal in STFT domain
            y_STFT = x_STFT(:,1,:)
            
        end
        
    end

    

%%
    %%% PROCESSING

    
    psi_STFT = calc_FD_GCC(y_STFT, dist_mat, phat_mode, weight_gcc_mode);
    
    
    if figure_mode
        fig = true_loc_idx;
    else
        fig = 0;
    end
    
    
    if prediction_3D_mode
        proj_dim = 3;
    else
        proj_dim = 2;
    end


    approxErr_dB = zeros(L, length(N_aux));
    locErr = zeros(L, length(N_aux));
    locErrRM = zeros(L, length(N_aux));
    DOAvec = DOAvec_list{1, end};

    
    [maxIdx_conv, index, maxIdx_conv_RM] = calc_SRPapprFast(psi_STFT, omega, T, N_mm, N_aux, Delta_t_list, DOAvec_list, fig);

    
    estim_DOAvec = DOAvec(maxIdx_conv,:, index);
    estim_DOAvec_RM = DOAvec(maxIdx_conv_RM,:, index);
    
    
    true_DOAvec_proj = zeros(1, 3);
    true_DOAvec_proj(1:proj_dim) = true_DOAvec(true_loc_idx,1:proj_dim);
    
    
    locErr(:) = rad2deg(acos(...
        (estim_DOAvec*transpose(true_DOAvec_proj))./(sqrt(sum(estim_DOAvec.^2, 2))*norm(true_DOAvec_proj))...
        ));
    
    locErrRM(:) = rad2deg(acos(...
        (estim_DOAvec_RM*transpose(true_DOAvec_proj))./(sqrt(sum(estim_DOAvec_RM.^2, 2))*norm(true_DOAvec_proj))...
        ));
    
    
    if graph_mode
        graph_err(true_loc_idx) = locErr;
        graph_errRM(true_loc_idx) = locErrRM;
    end


    %% SAVE

    res.approxErr_dB(true_loc_idx, :, :) = approxErr_dB;
    res.locErr(true_loc_idx, :, :) = locErr;
    res.locErrRM(true_loc_idx, :, :) = locErrRM;

    Tictoc(true_loc_idx) = toc;

end

%% STATISTICAL RESULTS


disp(['Mic array configuration: ' mic_config])
disp(['SNR: ' num2str(SNR)])
disp(['Num of targets: ' num2str(size(true_DOAvec, 1))])
disp(['mean time per target: ' num2str(mean(Tictoc))])
disp(['var time per target: ' num2str(var(Tictoc))])
disp(['mean error: ' num2str(mean(res.locErr))])
disp(['var error: ' num2str(var(res.locErr))])
disp(['mean RM error: ' num2str(mean(res.locErrRM))])
disp(['var RM error: ' num2str(var(res.locErrRM))])
disp('DONE.')


%%

if graph_mode
    
    if Attenuation_Mode
        atten_text = '_LPF_';
    else
        atten_text = '_';
    end
    
    if phat_mode
        phat_text = '_phated_';
    else
        phat_text = '_phatless_';
    end
    
    if weight_gcc_mode
        weight_text = 'weightled';
    else
        weight_text = 'weightless';
    end
    
    graph_file_name = [mic_config atten_text num2str(loc_num) '_' signal_type phat_text weight_text '.mat'];
    save(fullfile("data/graphs", graph_file_name), "graph_err", "graph_errRM", "graph_SNR", "graph_dist");
    
    error_graph(graph_file_name);
end
