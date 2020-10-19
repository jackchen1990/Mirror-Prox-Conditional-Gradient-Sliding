function x = FW_x(beta,x1,g1,eta,tau)

    maxiter = 1000;

    [d,h] = size(g1);


    x = x1;

    for k = 1:maxiter

        g = beta*(x-x1) + g1;
        
        if norm(g) < eta
            break;
        end
        
        if d > h
        g = g';
        end
        [u1,s1,~] = svds(g*g',1);
        v1 = g'*u1/sqrt(s1(1));
        y = -tau*u1*v1';
        if d > h
            y = y';
            g = g';
        end


        h = sum(sum((x-y).*g));
        if h <= eta
            break;
        end
        gamma = min(1, h / beta / norm(x-y,'fro')^2);
        x = (1-gamma)*x + gamma*y;

    end
    


end

