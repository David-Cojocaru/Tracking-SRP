function [ Psi_STFT ] = calc_FD_GCC(y_STFT, dist_mat, phat_mode, weight_mode)
% [ Psi_STFT ] = calc_FD_GCC(y_STFT)
% computes frequency domain GCCs.
%
% IN:
% y_STFT         microphone signal - frequency x frames x microphones
%
% OUT:
% % Psi_STFT     FD GCCs - frequency x frames x microphone pairs


[K,L,M] = size(y_STFT);
P = M*(M-1)/2;
Psi_STFT = zeros(K,L,P);
weight = 1 / min(dist_mat(find(dist_mat > 0)));

for k = 1:K
    
    for l = 1:L
        
        y = squeeze(y_STFT(k,l,:));

        psi=zeros(P,1); 
        p = 0;
        for mprime = 1:M
            for m = mprime+1:M
                p = p+1;
                psi(p) = y(m)*y(mprime)';
                
                if weight_mode
                    psi(p) = psi(p) * (weight * dist_mat(m, mprime)) ^ 10;
                end
            end
        end
       
        if phat_mode
            psi = psi./(abs(psi)+ 1e-9);
        end
        
        Psi_STFT(k,l,:) = psi;
        
    end
end

end