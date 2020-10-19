function [k,y] = FW_y1(beta,y1,g1,eta)

    maxiter = 100000;
    
    y = y1;
    
    for k = 1:maxiter

        g = beta*(y-y1) + g1;
        
        [~,ind] = min(g);
        
        s = zeros(size(g));
        s(ind) = 1;

        h = (y-s)' * g;
        
        %gamma = min(1, h / beta / norm(s-y)^2);
        

        
        if h <= eta
            break;
        end
        gamma = min(1, h / beta / norm(s-y)^2);
        
        y = y -gamma*(y-s);
        
        
        

    end
    
end
