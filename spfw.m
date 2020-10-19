function [time, g_fw] = spfw(A, L, tau, lambda, iter)

    g_fw = zeros(iter+1, 1);

    
    [n,d] = size(A);
    [~,h] = size(L);

    
    x = zeros(d,h);
    y = zeros(n,1);
    y(1) = 1;
    time = zeros(iter+1,1);
    
    [g_x,g_y,g_fw(1)] = compute_gap(A, L, x, y, tau, lambda);
        fprintf('gap = %f, g_x = %f, g_y = %f\n',g_fw(1), g_x, g_y);
    
        
    for t = 1:iter

        tic;
        
        gamma = 2 / (t+2);
        
        % compute g_fw

        grad_x = A' * (y .* (L - exp(A*x) ./ (sum(exp(A*x),2)))) / n;
        grad_y = - lambda * n * (n * y - ones(n, 1)) - sum(L .* (log(exp(A*x) ./ (sum(exp(A*x),2)))),2)/n;

        
        if d > h
            grad_x = grad_x';
        end
        [u1,s1,~] = svds(grad_x*grad_x',1);
        v1 = grad_x'*u1/sqrt(s1(1));
        

        s_x = -tau*u1*v1';
        if d > h
            s_x = s_x';
        end
        
        [~,ind] = min(-grad_y(:));
        [idx1,idx2] = ind2sub(size(grad_y),ind);
        s_y = zeros(size(grad_y));
        s_y(idx1,idx2) = 1;


        x = (1 - gamma) * x + gamma * s_x;
        y = (1 - gamma) * y + gamma * s_y;
        
        time(t+1) = toc;
        if t > 1
            time(t+1) = time(t+1) + time(t);
        end
        
        
       %% compute g_fw
       [g_x,g_y,g_fw(t+1)] = compute_gap(A, L, x, y, tau, lambda);
        fprintf('gap = %f, g_x = %f, g_y = %f\n',g_fw(t+1), g_x, g_y);
       
        
        
    end

end