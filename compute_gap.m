function [g_x,g_y,g_fw] = compute_gap(A, L, x, y, tau, lambda)

    [n,d] = size(A);
    [~,h] = size(L);
       
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
        grad_x = grad_x';
    end
    g_x = sum(sum((x-s_x).*grad_x));

    [~,ind] = min(-grad_y(:));
    [idx1,idx2] = ind2sub(size(grad_y),ind);
    s_y = zeros(size(grad_y));
    s_y(idx1,idx2) = 1;
    g_y = -(y-s_y)' * grad_y;

    g_fw = g_x + g_y;
end