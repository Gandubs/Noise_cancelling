w = load('ananth - external_noise.txt');     
y = load('ananth - noisy_speech.txt');        
s_clean = load('ananth - clean_speech.txt'); %will be used for SNR calculations

%Please dont choose M > batch_size
M = 5;        %gives best snr for given files (trials and error)     
batch_size = 1600; 
mu = 0.12;         

Wz = LMS(w, y, M, batch_size, mu);

v_hat = filter(Wz, 1, w);
s_hat = y - v_hat;  

compute_snr = @(clean, noisy) 10 * log10(sum(clean.^2) / sum((noisy - clean).^2));

snr_unfiltered = compute_snr(s_clean, y);  
snr_filtered  = compute_snr(s_clean, s_hat); 

fprintf('SNR before noise cancellation: %.2f dB\n', snr_unfiltered);
fprintf('SNR after noise cancellation: %.2f dB\n', snr_filtered);

disp('Playing the estimated clean speech...');
sound(s_hat, 44100);
