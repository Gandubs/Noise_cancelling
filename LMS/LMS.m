function [Wz, mse_values] = LMS(w, y, M, buffer_size, mu) % returned mse_values for analysis
% The LMS filter
% i/p:
%    w         - external noise vector (column vector) - i.e, w(n)
%    y         - final noisy speech (column vector), i.e., s(n) + v(n)
%    M         - order of the filter
%    buffer_size- pretty self explanatory
%    mu        - Learning rate
%
% o/p:
%    Wz        - filter coefficients (a "fitter" of sorts)
    N = length(w);
    num_batches = floor(N / buffer_size);
    Wz = zeros(M, 1);      
    mse_values = zeros(num_batches, 1); 
    for b = 1:num_batches
        start = (b-1) * buffer_size + 1;
        stop = start + buffer_size - 1;
        if stop > N
            break;
        end

        x_batch = w(start:stop); 
        d_batch = y(start:stop);    
        e_batch = zeros(buffer_size, 1);    
        for n = M:buffer_size 
            xvec = x_batch(n:-1:n-M+1);
            v_hat = Wz' * xvec;           
            e_batch(n) = d_batch(n) - v_hat;
            mod = xvec' * xvec;
            Wz = Wz + (mu/mod) * xvec * e_batch(n);
        end
        mse_values(b) = sum(e_batch.^2);
    end
end
