function band_w = filt_response(Fs, notch_freq, r)
%code only returns bandwidth if single notch freq is given, else 0.
    band_w = 0;
    
    %basically zoom in if we have one notch freq, else it will try to show
    %all notches in the best way
    if length(notch_freq) > 1
        freq_min = max(0, min(notch_freq) - 50);
        freq_max = min(Fs/2, max(notch_freq) + 50);
    else
        freq_min = max(0, notch_freq - 200);
        freq_max = notch_freq + 200;
    end
    
    N = 3000000; 
    w = linspace(0, Fs/2, N);
    b_cascade = 1;
    a_cascade = 1;

    %note that filter coeffs are updated using convolution
    for i = 1:length(notch_freq)
        w0 = 2 * pi * notch_freq(i) / Fs;
        b = [1, -2*cos(w0), 1];
        a = [1, -2*r*cos(w0), r*r];        
        b_cascade = conv(b_cascade, b);
        a_cascade = conv(a_cascade, a);
    end

    % freq response of final cascaded filter
    h_cascade = freqz(b_cascade, a_cascade, w, Fs);
    
    % Plot magnitude response in dB
    magnitude_db = 20 * log10(abs(h_cascade));
    figure('Position', [100, 100, 800, 600]);
    subplot(2, 1, 1);
    plot(w, magnitude_db, 'LineWidth', 1.5);
    grid on;
    
    if isscalar(notch_freq)
        title_str = sprintf('Notch Filter Magnitude Response (f₀ = %.3f Hz, r = %.5f)', notch_freq, r);      
        % We seach for the nearest sample(index) with mag<= -3, dies this for single
        % notch only btw
        idx = find(magnitude_db <= (max(magnitude_db)-3));
        if length(idx) >= 2
            f_low = w(idx(1));
            f_high = w(idx(end));
            band_w = f_high - f_low;
        end
    else
        title_str = sprintf('Cascaded Notch Filter Response (%d notches, r = %.5f)', length(notch_freq), r);
    end
    title(title_str);
    
    % Add reference lines
    hold on;
    ylimits = ylim;
    db_threshold = 20 * log10(1/sqrt(2)); % -3dB line
    plot([w(1), w(end)], [db_threshold, db_threshold], 'k--', 'LineWidth', 1);
    
    % add vertical lines at notch freqs
    colors = {'r--', 'g--', 'b--', 'm--', 'c--'};
    for i = 1:length(notch_freq)
        color_idx = mod(i-1, length(colors)) + 1;
        plot([notch_freq(i), notch_freq(i)], ylimits, colors{color_idx}, 'LineWidth', 1);
    end
    xlim([freq_min, freq_max]);
    
    % Plot phase response in degrees
    subplot(2, 1, 2);
    phase_deg = angle(h_cascade) * 180 / pi;
    plot(w, phase_deg, 'LineWidth', 1.5);
    grid on;
    title('Phase Response');
    xlabel('Frequency (Hz)');
    ylabel('Phase (degrees)');
    
    % Add notch frequency markers
    hold on;
    ylimits = ylim;
    for i = 1:length(notch_freq)
        color_idx = mod(i-1, length(colors)) + 1;
        plot([notch_freq(i), notch_freq(i)], ylimits, colors{color_idx}, 'LineWidth', 1);
    end
    xlim([freq_min, freq_max]);
end
