function [ x_TD, delay_diff_array ] = calc_INPUT_SIGNAL(x_TDmono, micPos, true_loc, c, fs, fs_new)

%% Mic Array Geometric Setup



% %%% MICROPHONE ARRAY
% % circular array, 10cm radius, six microphones 
% tmp = load('coord_mic_array.mat');
% % array center
% arrayCenterPos = tmp.arrayCenterPos;
% % microphone positions
% micPos = tmp.micPos;
% % number of microphones
% M = size(micPos,1);

% num_of_mics = 6;

% % anti-prism ensamble:
% mic_1 = [0.5, 0.5, 0];
% mic_2 = [-0.5, 0.5, 0];
% mic_3 = [-0.5, -0.5, 0];
% mic_4 = [0.5, -0.5, 0];
% mic_5 = [0.25, 0, 0.2];
% mic_6 = [0, 0.25, 0.2];
% mic_7 = [-0.25, 0, 0.2];
% mic_8 = [0, -0.25, 0.2];
% 
% mic_array = [mic_1; mic_2; mic_3; mic_4; mic_5; mic_6; mic_7; mic_8];

%% Source setup and 

% source = [100 0 0];
% source_array = [source; source; source; source; source; source; source; source];

% x_TDmono = resample(x_TDmono, fs_new, fs);

cord_diff_array = repmat(true_loc, size(micPos, 1), 1) - micPos;

% speed_of_sound = 340;

dist_diff_array = zeros(size(micPos, 1), 1);

for i = 1:size(micPos, 1)
    dist_diff_array(i) = norm(cord_diff_array(i, :));
end
% dist_diff_array = [norm(cord_diff_array(1, :)), norm(cord_diff_array(2, :)), norm(cord_diff_array(3, :)), norm(cord_diff_array(4, :)), norm(cord_diff_array(5, :)), norm(cord_diff_array(6, :)), norm(cord_diff_array(7, :)), norm(cord_diff_array(8, :))]
diff_array = dist_diff_array - min(dist_diff_array);
delay_diff_array = diff_array / c;
delay_diff_array_sample = floor(delay_diff_array * fs_new);

% x_TD = zeros(size(x_TDmono, 1), size(micPos, 1));

x_TDmono = resample(x_TDmono, fs_new, fs);

x_TD = delayseq(repmat(x_TDmono, 1, size(micPos, 1)), delay_diff_array_sample);





% for i = 1:size(micPos, 1)
%     x_TD
% end



%%

% Fs = 32000;
% [signal, Fs] = audioread('signal.wav');
% signal = signal(:, 1);
% 
% len = length(signal);

% winSize = round(Fs*len);
% overlap = round(winSize/2);
% fftSize = winSize;
% 
% spectrogram(signal, winSize, overlap, fftSize, Fs, 'yaxis');

% t = linspace(1, len, len);
% figure(1);
% plot(t, signal);
% 
% figure(2);
% spectrogram(signal, round(sqrt(len)), 'yaxis');






