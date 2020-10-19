function r = cgs(A, x, r, L, lambda, N, M, tau)
    
    [n,~] = size(A);
    %M = 5;
    for s = 1:N
        
        y = r;
        v = r;

        for k = 1:M
            
            alpha = 2/(k+1);
            z = (1-alpha) * y + alpha * v;
            grad_y = - lambda * n * (n * z - ones(n, 1)) - sum(L .* (log(exp(A*x) ./ (sum(exp(A*x),2)))),2)/n;
            
            beta = 2 * n / k;
            eta = 4 * n / k / M / 2^s;
            
            [it,v] = FW_y1(beta,v,-grad_y,eta);
            y = (1-alpha) * y + alpha * v;
            
            %if mod(k,4) == 0

                    %fprintf('y_gap = %f, it = %d, eta=%f\n',g_y, it,eta);
                    %fprintf(' it = %d, eta=%f\n',it,eta);
                
                %fprintf('gap = %f, g_x = %f, g_y = %f\n',g_fw, g_x, g_y);
            %end
            
        end
        r = y;
        %[g_x,g_y,g_fw] = compute_gap(A, L, x, y, tau, lambda);
        %fprintf('s = %d, fgap = %f, g_x = %f, g_y = %f\n',s, g_fw, g_x, g_y);
        
    end

end