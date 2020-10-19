function [x,y,u] = prox(A, L, tau, lambda, beta, eta, gamma_k, eps_t, x_k1, y_0,u, z)

[n,d] = size(A);
[~,h] = size(L);

x = x_k1;

eps_mp = 6*sqrt(2*eps_t);
eps_cgs = eps_t/32;
R = ceil(max(1,log2(1/eps_mp)));
T = ceil(max(1,log2(1/eps_cgs)));

for r1 = 1:R
    
    
    y = cgs(A, x, y_0, L, lambda, T, 5, tau);
    
    
    grad_x = A' * (y .* (L - exp(A*z) ./ (sum(exp(A*z),2)))) / n;

    u = FW_x(beta,u,grad_x,eta,tau);
    x = (1 - gamma_k) * x_k1 + gamma_k * u;

    
end


end