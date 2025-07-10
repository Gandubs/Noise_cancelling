function [proposed_metric_1, proposed_metric_2]  = partial_snr(partial_sup_op, s_clean, noisy, notch_freq)   
    fs = 44100;%Hz
    bw = 50;%Hz
    notched_op = partial_sup_op;
    notched_clean = s_clean;
    for i = 1:length(notch_freq)
        [b,a] = iirnotch(notch_freq(i)/(fs/2), bw/(fs/2));
        temp1 = filter(b,a,notched_op);
        notched_op = temp1;
        temp2 = filter(b,a,notched_clean);
        notched_clean = temp2;
    end
    proposed_metric_1 = 10 * log10(mean(notched_clean.^2) / mean((notched_op - notched_clean).^2));    
    %plotFFT(y, 44100, 'iirnotch applied to Partial Supression (for metric)');
    partial_sup_op= partial_sup_op(:);
    N = length(partial_sup_op); 
    X = fftshift(fft(noisy, N)); %fftshift centers the plot around 0
    f = (-N/2:N/2-1) * (fs / N);
    % magnitude in dB
    P = 20 * log10(abs(X));
    
    N_2 = length(partial_sup_op); 
    X_2 = fftshift(fft(partial_sup_op, N_2)); %fftshift centers the plot around 0
    f_2 = (-N_2/2:N_2/2-1) * (fs / N_2);
    % magnitude in dB
    P_2 = 20 * log10(abs(X_2));
    
    proposed_metric_2 = zeros(size(notch_freq));
    for i = 1:length(notch_freq)
        idx1 = find(f <= notch_freq(i), 1, 'last');
        max_mag1 = max(P(idx1 - round(2*N/fs):idx1 + round(2*N/fs)));
        idx2 = find(f_2 <= notch_freq(i), 1, 'last');
        max_mag2 = max(P_2(idx2 - round(2*N/fs):idx2 + round(2*N/fs)));
        proposed_metric_2(i) = max_mag1 - max_mag2;
    end
end