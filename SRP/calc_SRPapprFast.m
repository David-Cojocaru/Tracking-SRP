function [maxIdx_conv,  index, maxIdx_conv_RM] = calc_SRPapprFast(Psi_STFT,  omega, T, N_mm, N_aux, Delta_t_list, DOAvec_list, fig)
% [ SRP_stack, nT, xi_mm_samp ] = calc_SRPappr(Psi_STFT,  omega, T, N_mm, N_aux, Delta_t_i)
% computes SRP approximation. 
%
% IN:
% Psi_STFT        FD GCCs - frequency x frames x microphone pairs
% omega           frequency vector
% T               sampling period
% N_mm            number of samples inside TDOA interval
% N_aux           bumber of auxilary samples outside TDOA interval (if vector, results are computed for each element)
% Delta_t_i       TDOA vector
%
% OUT:
% SRP_stack       SRP map - frames x candidate locations

maxIdx_conv = 1;
maxIdx_conv_prev = maxIdx_conv;
size_prev = 0;

[K, L, P] = size(Psi_STFT);
xi_mm_samp = cell(P);

for i = 1:size(Delta_t_list, 2)
    
    Delta_t = Delta_t_list{1, i};
    DOAvec = DOAvec_list{1, i};

    index = (maxIdx_conv_prev - 1) * size_prev + maxIdx_conv;
    Delta_t_i = Delta_t(:, :, index);
    DOAvec_i = DOAvec(:, :, index);
    size_prev = size(DOAvec, 1);
    
%     disp('* compute SRP approximation...')

    [J, P] = size(Delta_t_i);

    % nT = cell(P);
    % xi_mm_samp = cell(P);

    SRP = zeros(J, 1);
    xi_mm_int = zeros(J, 1);

    for p = 1:P
        
            psi = Psi_STFT(:,1,p);

            % sample points
            N_half = N_mm(p) + N_aux;
            n = -N_half:1:N_half;
            nT = (n*T).';

        if i == 1
            xi_mm_samp{p} = real(exp(1j*nT*omega.')*Psi_STFT(:,1,p)); % (10); nDT{p}*omega.' is N x K matrix, psi is K x 1 vector

        end
        xi_mm_int(:) = sinc(Delta_t_i(:,p)/T - n(1:end))*xi_mm_samp{p}(1:end); % (14); Delta_t(:,p)/T - n is J x N matrix, xi_mm_samp is N x 1 vector
        SRP = SRP + 2*xi_mm_int; % (11)
    %     SRP = 0;
%         toc;

    end

    
    if fig ~= 0
        figure(fig);
        scatter3(DOAvec_i(:, 1), DOAvec_i(:, 2), DOAvec_i(:, 3), 70, transpose(SRP), 'filled');axis equal;
        cb = colorbar;
    end
    
    maxIdx_conv_prev = index;

    % localization error
    [~, maxIdx_conv] = max(SRP);
    maxIdx_conv_RM = 1;
    [estim_DOAvec_RM, maxIdx_conv_RM] = calc_RobustMax(SRP, DOAvec, 0.85, 7);

end

end
