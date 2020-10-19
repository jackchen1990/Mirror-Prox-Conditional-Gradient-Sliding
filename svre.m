function [time, g_fw] = svre(A, L, tau, lambda, eta_x, eta_y, iter, in_iter, batch_size)
    g_fw = zeros(iter+1, 1);
    
    [n,d] = size(A);
    [~,h] = size(L);
    
    x = rand(d,h);
    y = zeros(n,1);
    y(1) = 1;
    time = zeros(iter+1,1);
    

    for t = 1:iter
        
       %% update
        
       
        tic;
        
        x0 = x;
        y0 = y;
        grad_x = A' * (y0 .* (L - exp(A*x0) ./ (sum(exp(A*x0),2)))) / n;
        grad_y = - lambda * n * (n * y0 - ones(n, 1)) - sum(L .* (log(exp(A*x0) ./ (sum(exp(A*x0),2)))),2)/n;
        
        for k=1:in_iter
            x_old = x; y_old = y;
            
            idx = randperm(n,batch_size);
            
            gradi_x = A(idx,:)' * (y(idx) .* (L(idx,:) - exp(A(idx,:)*x) ./ (sum(exp(A(idx,:)*x),2)))) / batch_size;
            gradi_y = - lambda * n * (n * y - ones(n, 1));
            gradi_y(idx) = gradi_y(idx) - sum(L(idx,:) .* (log(exp(A(idx,:)*x) ./ (sum(exp(A(idx,:)*x),2)))),2) / batch_size;
            
            gradi_x0 = A(idx,:)' * (y0(idx) .* (L(idx,:) - exp(A(idx,:)*x0) ./ (sum(exp(A(idx,:)*x0),2)))) / batch_size;
            gradi_y0 = - lambda * n * (n * y0 - ones(n, 1));
            gradi_y0(idx) = gradi_y(idx) - sum(L(idx,:) .* (log(exp(A(idx,:)*x0) ./ (sum(exp(A(idx,:)*x0),2)))),2) / batch_size;
        
            x = x - eta_x * (gradi_x - gradi_x0 + grad_x);
            x = projection_x(x,tau);   
            y = y + eta_y * (gradi_y - gradi_y0 + grad_y);
            y = projection_y(y);
            
            idx = randperm(n,batch_size);
            
            gradi_x = A(idx,:)' * (y(idx) .* (L(idx,:) - exp(A(idx,:)*x) ./ (sum(exp(A(idx,:)*x),2)))) / batch_size;
            gradi_y = - lambda * n * (n * y - ones(n, 1));
            gradi_y(idx) = gradi_y(idx) - sum(L(idx,:) .* (log(exp(A(idx,:)*x) ./ (sum(exp(A(idx,:)*x),2)))),2) / batch_size;
            
            gradi_x0 = A(idx,:)' * (y0(idx) .* (L(idx,:) - exp(A(idx,:)*x0) ./ (sum(exp(A(idx,:)*x0),2)))) / batch_size;
            gradi_y0 = - lambda * n * (n * y0 - ones(n, 1));
            gradi_y0(idx) = gradi_y(idx) - sum(L(idx,:) .* (log(exp(A(idx,:)*x0) ./ (sum(exp(A(idx,:)*x0),2)))),2) / batch_size;
        
            x = x_old - eta_x * (gradi_x - gradi_x0 + grad_x);
            x = projection_x(x,tau);   
            y = y_old + eta_y * (gradi_y - gradi_y0 + grad_y);
            y = projection_y(y);

        end
        
        time(t+1) = toc;
        if t > 1
            time(t+1) = time(t+1) + time(t);
        end
        
        
       %% compute g_fw
       [g_x,g_y,g_fw(t+1)] = compute_gap(A, L, x, y, tau, lambda);
        fprintf('gap = %f, g_x = %f, g_y = %f\n',g_fw(t+1), g_x, g_y);
        
       
     end
end