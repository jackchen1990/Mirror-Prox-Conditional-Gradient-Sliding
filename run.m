clc;
clear;

load('dataset/news20');

[n, d] = size(A);
if min(b) == 0
    b = b + 1;
end

h = max(b);
L = zeros(n,h);
for i = 1:n
    L(i,b(i)) = 1;
end

tau = 100;
lambda = 1/n;

batch = 100;
in_iter = ceil(n/batch);

fprintf('spfw\n');
iter = 600;
[time_spfw,g_fw_spfw] = spfw(A, L, tau, lambda, iter);


fprintf('svre\n');
iter = 5;
eta_x = 1e-1;
eta_y = 1e-1;
[time_svre, g_fw_svre] = svre(A, L, tau, lambda, eta_x, eta_y, iter, in_iter, batch);

iter = 20;

fprintf('mpscgs\n');
[time_mpscgs,g_fw_mpscgs] = mpscgs(A, L, tau, lambda, iter);
fprintf('mpcgs\n');
[time_mpcgs,g_fw_mpcgs] = mpcgs(A, L, tau, lambda, iter);

save('news20.mat','time_svre','g_fw_svre','time_spfw','g_fw_spfw','time_mpscgs','g_fw_mpscgs','time_mpcgs','g_fw_mpcgs');

