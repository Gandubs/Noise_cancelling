    function s_hat = Full_Supp(x, d, lambda, delta, M)
        s_hat = zeros(length(x),1);
        n = 1;
        while n <= M-1
            s_hat(n) = x(n);
            n = n + 1;
        end
        P = (1/delta)*eye(M);
        Wz = zeros(M,1);

        while n <= length(x)
            xvec = x(n:-1:n-M+1); %external noise
            z = P * xvec;  % p(n-1)xvec(n)
            g = (1/(lambda + xvec' * z))*z; 
            alpha = d(n) - Wz' * xvec; % priori error
            Wz = Wz + alpha*g; % weight update eqn
            P = (1/lambda)*(P - g * xvec' * P); % updation of inverse correln matrix 
            s_hat(n) = d(n) - (Wz' * xvec); % update 'predicted' signal
            n = n + 1;
        end

    end