function x = projection_x(x,tau)
    [u,s,v] = svd(x,'econ');
    s = diag(s);
    ss = 0;
    for i = 1:length(s)
        ss = ss + s(i);
        if ss > tau
            s(i) = s(i) - (ss - tau);
            s(i+1:end) = 0;
            break;
        end
    end
    x = u * diag(s) * v';

end