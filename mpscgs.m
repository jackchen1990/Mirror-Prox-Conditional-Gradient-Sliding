function [time, g_fw] = mpscgs(A, L, tau, lambda, iter)

    g_fw = zeros(iter+1, 1);
    time = zeros(iter+1,1);
    
    [n,d] = size(A);
    [~,h] = size(L);

    sum_y = 0;
    x = zeros(d,h);
    y = zeros(n,1);
    y(1) = 1;
    y_bar = y;
    v = x;
    
    for k = 1:iter
        
        tic;
        
        gamma_k = 3 / (k+2);
        beta = 6*n / (k+1);
        eta = n / 24 / k / (k+1);
        z = (1-gamma_k) * x + gamma_k * v;
        
        eps_t = 1/k/(k+1)/(k+2);
        
        B = min(n,k^3);
        idx = randperm(n,B);
        
        [x,y,v] = sprox(A, L, tau, lambda, beta, eta,gamma_k, eps_t, x, y,v, z,idx,B);
        sum_y = sum_y + y*k*(k+1);

        y_bar = 3*sum_y / (k+1) / k / (k+2);
        
        
        time(k+1) = toc;
        if k > 1
            time(k+1) = time(k+1) + time(k);
        end
        
        
       %% compute g_fw
       [g_x,g_y,g_fw(k+1)] = compute_gap(A, L, x, y_bar, tau, lambda);
        fprintf('gap = %f, g_x = %f, g_y = %f\n',g_fw(k+1), g_x, g_y);
        
        
        
    end
end