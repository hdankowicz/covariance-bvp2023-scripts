function pts = findPoints(data, vpt)
% FINDPOINTS     Locate point on torus by its parameterization

npt = size(vpt,1);
rho = data.rho;
for i=1:npt
  if vpt(i,2)>=1 || vpt(i,2)<0
    itau = floor(vpt(i,2));
    vpt(i,2) = vpt(i,2)-itau;
    vpt(i,1) = vpt(i,1)+rho*itau;
  end
end

F1d = data.F1d;
N   = data.N;
NTST = data.NTST;
NCOL = data.NCOL;

jbp = floor(interp1(0:1/NTST:1,1:NTST+1,vpt(:,2),'linear'));
rows = repmat((jbp(:)-1)'*(NCOL+1), [(NCOL+1),1]);
rows = rows+repmat((1:(NCOL+1))',[1, npt]);
cols = repmat(1:npt,[NCOL+1,1]);
tc   = 2*NTST*vpt(:,2)+1-2*jbp;
tm   = linspace(-1,1,NCOL+1);
Lc   = coll_L(tm, tc);
Lbp  = full(sparse(rows,cols,Lc',(NCOL+1)*NTST,npt));
Wmap = kron(Lbp',eye(4));
Xtau = Wmap*data.ptmesh;

pts = zeros(4,npt);
for i=1:npt
  cs = cos(2*pi*(1:N).*vpt(i,1));
  ss = sin(2*pi*(1:N).*vpt(i,1));
  B = [1 reshape([cs; ss], [1 2*N])];
  xtau = reshape(Xtau(4*(i-1)+(1:4),:), [4*(2*N+1) 1]);
  pts(:,i) = reshape(kron(B*F1d,eye(4))*xtau,[4 1]);
end

end

function A = coll_L(ts, tz)

q = numel(ts);
p = numel(tz);

zi = repmat(reshape(tz, [p 1 1]), [1 q q]);
sj = repmat(reshape(ts, [1 q 1]), [p 1 q]);
sk = repmat(reshape(ts, [1 1 q]), [p q 1]);

t1 = zi-sk;
t2 = sj-sk;
idx = find(abs(t2)<=eps);
t1(idx) = 1;
t2(idx) = 1;

A = prod(t1./t2, 3);

end
