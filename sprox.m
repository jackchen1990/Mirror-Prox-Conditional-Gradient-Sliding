function [x,y,u] = sprox(A, L, tau, lambda, beta, eta, gamma_k, eps_t, x_k1, y_0,u, z,idx,batch_size)

[n,d] = size(A);
[~,h] = size(L);

x = x_k1;

eps_mp = 6*sqrt(2*eps_t);
eps_cgs = eps_t/64;
R = ceil(max(1,log2(1/eps_mp)));
T = ceil(max(1,log2(1/eps_cgs)));

for r1 = 1:R
    
    
    y = storc_dro(A, x, y_0, L, lambda, T, 5, 100,tau);
    gradi_x = A(idx,:)' * (y(idx) .* (L(idx,:) - exp(A(idx,:)*z) ./ (sum(exp(A(idx,:)*z),2)))) / batch_size;
    
    
    u = FW_x(beta,u,gradi_x,eta,tau);
    x = (1 - gamma_k) * x_k1 + gamma_k * u;


end


end