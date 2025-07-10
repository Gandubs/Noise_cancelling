function plotFFT(x, Fs, plotTitle)
    x = x(:);
    N = length(x); 
    
    X = fftshift(fft(x, N)); 
    f = (-N/2:N/2-1) * (Fs / N);
    P = 20 * log10(abs(X));

    figure;
    plot(f, P, 'LineWidth', 1.2);
    grid on;
    xlim([-Fs/2 Fs/2]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title(plotTitle);
end
