function s_hat = Partial_Supp(x, d, lambda, delta, M, r, notch_freq)
        s_hat = zeros(length(x),1);
        xvec_filt = zeros(M,1);
        Fs = 44100; 
        %filter coeffs : b = [b(1) b(2) b(3)] , a = [a(1) a(2) a(3)]
        %using MATLAB indexing btw
        n = 1;
        buffer_x_matrix = zeros(2, length(notch_freq));
        buffer_y_matrix = zeros(2, length(notch_freq));
        %Processes only first M samples in filt loop
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
            n = n + 1;
        end

        P = (1/delta)*eye(M);
        Wz = zeros(M,1);  

        while n <= length(x)
            %depending on how many notch filters are chosen, those many
            %times the input signal is looped through, and the resulting
            %filtered signal is passed to rls
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
            xvec = flip(xvec_filt); %same equations as used in full_supp
            z = P * xvec;
            g = (1/(lambda + xvec' * z))*z;
            alpha = d(n) - Wz' * xvec;
            Wz = Wz + alpha*g;
            P = (1/lambda)*(P - g * xvec' * P);
            s_hat(n) = d(n) - (Wz' * xvec);
            n = n + 1;
        end

end
