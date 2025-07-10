%used this file to examine characterisitcs of noise as a possible
%explanation for weird convergence behaviour, however no conclusion was drawn
w = load('K Rahul - external_noise.txt'); 
fs = 44100; 

% Compute fft
N = length(w);
W = abs(fft(w));
f = (0:N-1) * fs / N; 
% Plot fft
figure;
plot(f(1:floor(N/2)), 20*log10(W(1:floor(N/2))), 'r'); 
title('FFT of external noise');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;

y = load('K Rahul - noisy_speech.txt'); 
% Compute fft
Ny = length(y);
Y = abs(fft(y));
f_y = (0:Ny-1) * fs / Ny; 
% Plot fft
figure;
plot(f_y(1:floor(Ny/2)), 20*log10(Y(1:floor(Ny/2))), 'm'); 
title('fft of noisy speech');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;

s_clean = load('K Rahul - clean_speech.txt'); % True clean speech
%compte fft
Nc = length(s_clean);
S_clean = abs(fft(s_clean));
f_c = (0:Nc-1) * fs / Nc; 
%-plot fft
figure;
plot(f_c(1:floor(Nc/2)), 20*log10(S_clean(1:floor(Nc/2))), 'g'); 
title('fft of clean speech');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;


%s_hat = load('K Rahul - estimated_clean_speech.txt'); 
%compute fft
%Ns = length(s_hat);
%S_hat = abs(fft(s_hat));
%f_s = (0:Ns-1) * fs / Ns; 
% plot fft
%figure;
%plot(f_s(1:floor(Ns/2)), 20*log10(S_hat(1:floor(Ns/2))), 'b'); 
%title('fft of est clean speech');
%xlabel('Frequency (Hz)');
%ylabel('Magnitude (dB)');
%grid on;
