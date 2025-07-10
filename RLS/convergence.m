function convergence(w, y, lambda, delta, M, lambda_p, delta_p, M_p, notch_freq, r)
    % The same RLS code but now modified to return errors as well for
    % plotting
    [~, errors_full] = Full_Supp_with_errors(w, y, lambda, delta, M);
    [~, errors_partial] = Partial_Supp_with_errors(w, y, lambda_p, delta_p, M_p, notch_freq, r);
    
    %considering squared errors
    errors_full_squared = errors_full.^2;
    errors_partial_squared = errors_partial.^2;
    
    % Apply downsampling
    % keep every 300th sample, did this because convergence plot was looking too cluttered  
    samples = 1:length(errors_full);
    downsampled_indices = 1:300:length(errors_full);
    downsampled_samples = samples(downsampled_indices);
    downsampled_errors_full = errors_full_squared(downsampled_indices);
    downsampled_errors_partial = errors_partial_squared(downsampled_indices);
    figure('Position', [100, 100, 800, 500]);
    
    % Plot squared errors
    plot(downsampled_samples, downsampled_errors_full, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(downsampled_samples, downsampled_errors_partial, 'g-', 'LineWidth', 1.5);
    
    grid on;
    xlabel('Sample Number');
    ylabel('Squared Error e(n)^2');
    title('Error Convergence');
    legend('Full Suppression', 'Partial Suppression');
    set(gca, 'YScale', 'log');
    
    % algo that returns with error values
    function [s_hat, errors] = Full_Supp_with_errors(x, d, lambda, delta, M)
        s_hat = zeros(length(x),1);
        errors = zeros(length(x),1);  
        n = 1;
        while n <= M-1
            s_hat(n) = x(n);
            errors(n) = 0;  
            n = n + 1;
        end
        
        P = (1/delta)*eye(M);
        Wz = zeros(M,1);

        while n <= length(x)
            xvec = x(n:-1:n-M+1);
            z = P * xvec;
            g = (1/(lambda + xvec' * z))*z;
            err = d(n) - Wz' * xvec;  
            errors(n) = err; %error computation, the only change done from
                             %prev code.
            Wz = Wz + err*g;
            P = (1/lambda)*(P - g * xvec' * P);
            s_hat(n) = d(n) - (Wz' * xvec);
            if isnan(s_hat(n))
                fprintf("At n = %d, it becomes NaN", n);
            end
            n = n + 1;
        end
    end

    function [s_hat, errors] = Partial_Supp_with_errors(x, d, lambda, delta, M, notch_freq, r)
        s_hat = zeros(length(x),1);
        errors = zeros(length(x),1);  
        xvec_filt = zeros(M,1);
        Fs = 44100; 
        n = 1;
        buffer_x_matrix = zeros(2, length(notch_freq));
        buffer_y_matrix = zeros(2, length(notch_freq));
        while n <= M
            for i = 1:length(notch_freq)
                w0 = 2*pi*notch_freq(i)/Fs;
                b = [1, -2*cos(w0), 1];
                a = [1, -2*r*cos(w0), r*r];
                if i==1
                    x_rn = x(n);
                else
                    x_rn = xvec_filt(end);
                end
                x_filt_rn = b(1)*x_rn + b(2)*buffer_x_matrix(2,i) + b(3)*buffer_x_matrix(1,i)...
                - a(2)*buffer_y_matrix(2,i) - a(3)*buffer_y_matrix(1,i);
                if i==1
                    xvec_filt = [xvec_filt(2:M);x_filt_rn];
                else
                    xvec_filt(end) = x_filt_rn;
                end
                buffer_x_matrix(:,i) = [buffer_x_matrix(2,i);x_rn];
                buffer_y_matrix(:,i) = [buffer_y_matrix(2,i); x_filt_rn];
            end
            errors(n) = 0;  
            n = n + 1;
        end

        P = (1/delta)*eye(M);
        Wz = zeros(M,1);  

        while n <= length(x)
            for i = 1:length(notch_freq)
                w0 = 2*pi*notch_freq(i)/Fs;
                b = [1, -2*cos(w0), 1];
                a = [1, -2*r*cos(w0), r*r];
                if i==1
                    x_rn = x(n);
                else
                    x_rn = xvec_filt(end);
                end
                x_filt_rn = b(1)*x_rn + b(2)*buffer_x_matrix(2,i) + b(3)*buffer_x_matrix(1,i)...
                - a(2)*buffer_y_matrix(2,i) - a(3)*buffer_y_matrix(1,i);
                if i==1
                    xvec_filt = [xvec_filt(2:M);x_filt_rn];
                else
                    xvec_filt(end) = x_filt_rn;
                end
                buffer_x_matrix(:,i) = [buffer_x_matrix(2,i);x_rn];
                buffer_y_matrix(:,i) = [buffer_y_matrix(2,i); x_filt_rn];
            end
            xvec = flip(xvec_filt);
            z = P * xvec;
            g = (1/(lambda + xvec' * z))*z;
            err = d(n) - Wz' * xvec;
            errors(n) = err;
            Wz = Wz + err*g;
            P = (1/lambda)*(P - g * xvec' * P);
            s_hat(n) = d(n) - (Wz' * xvec);
            if isnan(s_hat(n))
                fprintf("At n = %d, it becomes NaN\n", n);
                break;
            end
            n = n + 1;
        end
    end
end