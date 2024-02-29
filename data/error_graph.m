function [] = error_graph(graph_file_name)

tmp = load(fullfile("data/graphs", graph_file_name));
% circ_rad = tmp.circ_rad;
% circ_err_mean = tmp.circ_err_mean
% circ_err = tmp.circ_err;
% circ_snr = tmp.SNR_dyn_mat;
graph_err = tmp.graph_err;
graph_SNR = tmp.graph_SNR;
graph_dist = tmp.graph_dist;

% figure(1)
% hold on
% for i = 1:size(circ_rad, 2)
%     p = plot(circ_snr(i, :), circ_err(i, :), '.', 'DisplayName', ['circle with radius of ' num2str(circ_rad(i))]);
%     [SNR_sorted, ind_sorted] = sort(circ_snr(i, :), 'descend');
%     err_sorted = circ_err(i, ind_sorted);
%     N = 50;
% %     mat = [circ_snr(i, :); circ_err(i, :)];
%     arr = movmean(err_sorted, N);
%     plot(SNR_sorted, arr, 'LineWidth', 2, 'Color', p.Color, 'DisplayName', ['circle with radius of ' num2str(circ_rad(i)) ' Moving Average']);
% end
% legend()
% xlabel('SNR (dB)')
% ylabel('Error (Deg)')
% title('UAV')
% hold off

figure(1)
hold on

p = plot(graph_SNR, graph_err, '.', 'DisplayName', ['platform']);
[SNR_sorted, ind_sorted] = sort(graph_SNR, 'descend');
err_sorted = graph_err(ind_sorted);
N = 50;
%     mat = [circ_snr(i, :); circ_err(i, :)];
arr = movmean(err_sorted, N);
plot(SNR_sorted, arr, 'LineWidth', 2, 'Color', p.Color, 'DisplayName', ['platform Moving Average']);

legend()
xlabel('SNR (dB)')
ylabel('Error (Deg)')
title('UAV')
hold off

% tmp = load(fullfile("data/graphs", 'circ_pairs_LPF_100_sig'));
% circ_rad = tmp.circ_rad;
% circ_err_mean = tmp.circ_err_mean
% circ_err = tmp.circ_err;
% circ_snr = tmp.SNR_dyn_mat;
% 
% figure(2)
% hold on
% for i = 1:size(circ_rad, 2)
%     plot(circ_snr(i, :), circ_err(i, :), '.', 'DisplayName', ['circle with radius of ' num2str(circ_rad(i))]);
% end
% legend()
% xlabel('SNR (dB)')
% ylabel('Error (Deg)')
% title('Explosion')
% hold off
