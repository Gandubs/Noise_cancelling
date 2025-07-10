% Load audio data
w       = load('ananth - external_noise.txt');
y       = load('ananth - noisy_speech.txt');
s_clean = load('ananth - clean_speech.txt');
fs = 44100;

%Full_suppression parameters
M       = 5;
lambda  = 0.999999;
delta   = 0.001;

% Partial_suppression parameters
M_p     = 5;
lambda_p = 0.999999;
delta_p = 0.001;

% Filter_Parameters
notch_freq = 999.9;       
r = 0.999;  % little phase and good proposal 1 metric, decent snr
%r=0.95;
%r = 0.999887; %very less phase but bad proposal 1 metric for notch at 999.9Hz

% Compute SNR before and after full suppression
s_hat = Full_Supp(w, y, lambda, delta, M);
snr_unfiltered = 10 * log10(mean(s_clean.^2) / mean((y - s_clean).^2));
snr_filtered   = 10 * log10(mean(s_clean.^2) / mean((s_hat - s_clean).^2));

fprintf('\nFull Suppression:\n');
fprintf('SNR before cancellation: %.2f dB\n', snr_unfiltered);
fprintf('SNR after  cancellation: %.2f dB\n', snr_filtered);
audiowrite('Full_Suppression.wav', s_hat, fs);

%Compute SNR after partial suppression
s_hat2 = Partial_Supp(w, y, lambda_p, delta_p, M_p, r, notch_freq);
snr_partial = 10 * log10(mean(s_clean.^2) / mean((s_hat2 - s_clean).^2));
fprintf('\nPartial Suppression:\n');
fprintf('Applied Notch frequencies: ');
fprintf('%.2f Hz ', notch_freq);
fprintf('\nSNR after cancellation: %.2f dB\n', snr_partial);
audiowrite('Partial_Suppression.wav', s_hat2, fs);

%Save Noisy Signal
audiowrite('Noisy_Speech.wav', y, fs);

% Plot FFTs
plotFFT(y, fs, 'Noisy Speech');
plotFFT(s_hat, fs, 'Full Suppression');
plotFFT(s_hat2, fs, 'Partial Suppression');

% Obtain Convergence plots
convergence(w, y, lambda, delta, M, lambda_p, delta_p, M_p, notch_freq, r);

% Gives the cascaded filter response, and return the bandwidth if
%only 1 notch frequency is input, else returns 0.
band_w = filt_response(fs, notch_freq, r);
fprintf('3 dB bandwidth: %.4f Hz\n', band_w);

%plotFFT(w,fs,'External Noise')

%Our proposed metrics for partial suppression
[proposed_metric_1, proposed_metric_2] = partial_snr(s_hat2, s_clean, y, notch_freq);
fprintf('Proposal 1: %.4f dB\n', proposed_metric_1);
fprintf('Proposal 2: %.4f dB\n', proposed_metric_2);


