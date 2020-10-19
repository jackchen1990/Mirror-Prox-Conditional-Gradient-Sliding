function y = projection_y(y)
    y = max(y,0);
    y = y / norm(y,1);
end