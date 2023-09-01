function [data, y] = projection(~, data, u)
% PROJECTION   Zero problem for identifying coordinates of torus projection

pt  = u(1:4);
vpt = u(5:6)';

data.ptmesh = data.Xbp_rep;
xpt = findPoints(data, vpt);

data.ptmesh = data.Lbp_rep1;
l1pt = findPoints(data, vpt);

data.ptmesh = data.Lbp_rep2;
l2pt = findPoints(data, vpt);

yT   = sum(l1pt.*(pt-xpt),1);
yrho = sum(l2pt.*(pt-xpt),1);

y = [yT(:);yrho(:)];

end