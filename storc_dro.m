function r = storc_dro(A, x, r, L, lambda, N, M, B,tau)
    
    [n,~] = size(A);

    for t = 1:N
        y = r;
        v = r;
        y_0 = y;
        
        grad_y = - lambda * n * (n * y - ones(n, 1)) - sum(L .* (log(exp(A*x) ./ (sum(exp(A*x),2)))),2)/n;
        
        for k = 1:M
            
            alpha = 2/(k+1);
            z = (1-alpha) * y + alpha * v;
            
            idx = randperm(n,B);
            
            gradi_z = - lambda * n * (n * z - ones(n, 1));
            gradi_z(idx) = gradi_z(idx) - sum(L(idx,:) .* (log(exp(A(idx,:)*x) ./ (sum(exp(A(idx,:)*x),2)))),2)/B;
            gradi_y = - lambda * n * (n * y_0 - ones(n, 1));
            gradi_y(idx) = gradi_y(idx) - sum(L(idx,:) .* (log(exp(A(idx,:)*x) ./ (sum(exp(A(idx,:)*x),2)))),2)/B;
            
            g = gradi_z - gradi_y + grad_y;
            
            
            
            beta = 2*n / k;
            eta = n / k / M / 2^t;
            
            [it,v] = FW_y(beta,v,-g,eta);
            y = (1-alpha) * y + alpha * v;
            
        end
        r = y;
        
    end

end