function [data, y] = interpolant(~, data, u)
% INTERPOLANT    Straight-line interpolant between data.x1 and data.x2
pt  = u(1:4);
ell = u(5);

y = pt-data.x1-ell*(data.x2-data.x1);

end
