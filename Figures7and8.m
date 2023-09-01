%% Analysis of noise-perturbed quasiperiodic invariant torus
%
% Script generating Figures 7 and 8 in Z. Ahsan, H. Dankowicz, C. Kuehn,
% "Adjoint-Based Projections for Uncertainty Quantification near
% Stochastically Perturbed Limit Cycles and Tori"

%% Process pre-computed data from COCO
% Run Figure6.m script first if error is thrown.

data = struct();

chart = coco_read_solution('rho1', 'two_tori_with_adjoints', 1, 'chart');
rho = chart.x;
data.rho  = rho;

[sol, bdata]= bvp_read_solution('torus1', 'two_tori_with_adjoints', 1);
nsegs = bdata.nsegs;
F1d   = bdata.bc_data.F1d;
N     = (nsegs-1)/2;
data.F1d  = F1d;
data.N    = N;

[~, cdata] = coll_read_solution('torus1.bvp.seg1', ...
  'two_tori_with_adjoints', 1);
NTST = cdata.coll_seg.maps.NTST;
NCOL = cdata.coll_seg.int.NCOL;
tbp  = cdata.coll_seg.mesh.tbp;
Xbp_rep = zeros(NTST*(NCOL+1)*4,nsegs);
for i=1:nsegs
  xbp_rep = interp1(sol{i}.tbp/sol{i}.tbp(end), sol{i}.xbp, tbp)';
  Xbp_rep(:,i) = xbp_rep(:);
end
data.NTST = NTST;
data.NCOL = NCOL;

Lbp_rep1 = zeros(NTST*(NCOL+1)*4,nsegs);
for i=1:nsegs
  fid = sprintf('torus1.bvp.seg%d',i);
  sol = coll_read_adjoint(fid, 'two_tori_with_adjoints', 1);
  lbp_rep = interp1(sol.tbp/sol.tbp(end), sol.l0', tbp)';
  Lbp_rep1(:,i) = lbp_rep(:);
end

Lbp_rep2 = zeros(NTST*(NCOL+1)*4,nsegs);
for i=1:nsegs
  fid = sprintf('torus2.bvp.seg%d',i);
  sol = coll_read_adjoint(fid, 'two_tori_with_adjoints', 1);
  lbp_rep = interp1(sol.tbp/sol.tbp(end), sol.l0', tbp)';
  Lbp_rep2(:,i) = lbp_rep(:);
end

covar = coco_bd_col('two_tori_with_adjoints', 'covariance');
Cbp_rep = zeros(NTST*(NCOL+1)*16,nsegs);
for i=1:nsegs
  Cbp_rep(:,i) = reshape(covar(:,:,:,i),NTST*(NCOL+1)*16,1);
end

data.Xbp_rep = Xbp_rep;
data.Lbp_rep1 = Lbp_rep1;
data.Lbp_rep2 = Lbp_rep2;
data.Cbp_rep = Cbp_rep;

%% run simulation
T  = 5000;
dt = 10^(-4);
N = 11600;

x1 = 0.1121;
x2 = 1.4771;
x3 = -0.2897;
x4 = 3.7509;

Xem = zeros(20000000,4);

% Simulation setup
tic
j = 1;
a = 1;
b = 1;
in = 1;
int = zeros(N,2);
k=1;
sigma   = 0.1;
epsilon = 0.5;
beta    = 0.5;
delta   = 1.9422;
while k<=N
  Fnoise = sqrt(dt)*randn(1,2);
  
  f = [x2; -epsilon*(x1^2-1)*x2-x1+beta*(x3-x1);
    x4; -epsilon*(x3^2-1)*x4-(1+delta)*x3+beta*(x1-x3)];
  
  x1 = x1+dt*f(1);
  x2 = x2+dt*f(2)+sigma*Fnoise(1);
  x3 = x3+dt*f(3);
  x4 = x4+dt*f(4)+sigma*Fnoise(2);
  if x1<0.4 && x1>-0.6 && x3<1-4.5*(x1+.2) && x3>1-4.5*(x1+.4)
    if ~in
      in = 1;
      int(k,:) = [a b];
      k = k+1;
      a = j;
    end
    b = j;
    Xem(j,:) = [x1 x2 x3 x4];
    j = j+1;
  else
    in = 0;
  end
end
toc
Xem = Xem(1:j-1,:);

save("stochastic_quasiperiodic_ex2_11600_05_raw.mat", "int", "Xem")

%% Find intersections

load('stochastic_quasiperiodic_ex2_11600_05_raw.mat');

as = [];
for k=2:2:size(int,1)
  k
  a = int(k,1); b = int(k,2);
  while b-a>1
    ind = floor((b+a)/2);
    dist = pdist2(Xem(ind,:),reshape(Xbp_rep,[4 numel(Xbp_rep)/4])');
    [val,idx_closest] = min(dist,[],2);
    Phi0 = (ceil(idx_closest/length(tbp))-1)/nsegs;
    tau0_idx = idx_closest-floor(idx_closest/length(tbp))*length(tbp);
    tau0_idx(tau0_idx==0) = length(tbp);
    tau0 = tbp(tau0_idx);
    prob = coco_prob;
    prob = coco_add_func(prob,'test', @projection, data, 'zero', ...
      'u0', [Xem(ind,:)'; [Phi0; tau0]]);
    prob = coco_add_pars(prob, 'pars', 1:4, {'x1','x2','x3','x4'});
    prob = coco_set(prob, 'cont', 'LogLevel', 0);
    coco(prob, 'test', [], 0)
    chart = coco_read_solution('test','test',1,'chart');
    if chart.x(6)<0.5
      a = ind;
    else
      b = ind;
    end
  end
  as = [as; a chart.x(5:6)'];
end

%% Compute points of intersection

ptsec = zeros(6,size(as,1));
for i=1:size(as,1)
  i
  prob = coco_prob;
  prob = coco_add_func(prob,'test', @projection, data, 'zero', ...
    'u0', [Xem(as(i,1),:)'; as(i,2:3)']);
  prob = coco_add_pars(prob, 'pars', 6, 'tau');
  prob = coco_set_parival(prob, 'tau', 0.5);
  idata = [];
  idata.x1 = Xem(as(i,1),:)';
  idata.x2 = Xem(as(i,1)+1,:)';
  prob = coco_add_func(prob, 'interp', @interpolant, idata, 'zero', ...
    'uidx', 1:4, 'u0', 0);
  prob = coco_set(prob, 'cont', 'LogLevel', 0);
  coco(prob, 'cross', [], 0)
  chart = coco_read_solution('test','cross',1,'chart');
  ptsec(:,i) = chart.x;
end

data.ptmesh = data.Xbp_rep;
xpt = findPoints(data, ptsec(5:6,:)');

save("stochastic_quasiperiodic_ex2_11600_05_processed.mat", "xpt", "ptsec")

%% load data
load("stochastic_quasiperiodic_ex2_11600_05_processed.mat")

%% Figure 7

figure
box on
grid on
hold on
thm = struct();
thm.sol.RO = {'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Marker', '.', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 12};
coco_plot_sol(thm, 'two_tori', 4, 'torus1', 'x', 'x', 'x')
plot3(xpt(1,:),xpt(2,:),xpt(3,:),'ro')
plot3(ptsec(1,:),ptsec(2,:),ptsec(3,:),'b.')
axis tight
view(-100,15)
set(gca, 'LineWidth', 2, 'Fontsize', 14)
xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 20)
ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 20)
zlabel('$x_3$', 'Interpreter', 'Latex', 'FontSize', 20)
hold off

%% Figure 8

data.ptmesh = data.Cbp_rep;
Cpt = reshape(findCovars(data, ptsec(5:6,:)'), [4,4,5800]);

[vs,es] = eigs(Cpt(:,:,1));
[es,js] = sort(diag(es));
ref3 = vs(:,js(3));
ref4 = vs(:,js(4));
dps = zeros(size(ptsec,2),2);
for i=1:size(ptsec,2)
  [vs,es] = eigs(Cpt(:,:,i));
  [es,js] = sort(diag(es));
  dx      = ptsec(1:4,i)-xpt(:,i);
  if ref3'*vs(:,js(3))<0
    vs(:,js(3)) = -vs(:,js(3));
  end
  if ref4'*vs(:,js(4))<0
    vs(:,js(4)) = -vs(:,js(4));
  end
  dps(i,:) = [dx'*vs(:,js(3)), dx'*vs(:,js(4))];
end

figure(4)
clf
hold on
plot(mod(ptsec(5,:),1), dps(:,1), 'r.')
grid on
box on
axis tight

figure(5)
clf
hold on
plot(mod(ptsec(5,:),1), dps(:,2), 'r.')
grid on
box on
axis tight

Crv = reshape(findCovars(data, [0:.01:1; 0.5*ones(1,101)]'), [4,4,101]);
vals = zeros(101,2);
for i=1:101
  eig = sort(eigs(Crv(:,:,i)));
  vals(i,:) = 0.1*sqrt(eig(3:4));
end
figure(4)
for radius=[-2,-1,1,2]
  plot(0:.01:1, radius*vals(:,1), 'k', 'LineWidth', 2);
end
plot(0:.01:1, zeros(101,1), 'k--', 'LineWidth', 2)
figure(5)
for radius=[-2,-1,1,2]
  plot(0:.01:1, radius*vals(:,2), 'k', 'LineWidth', 2);
end
plot(0:.01:1, zeros(101,1), 'k--', 'LineWidth', 2)

for m=0:1/29:1
  subY   = ptsec(:,abs(ptsec(5,:)-m)<1/58 | abs(ptsec(5,:)-m-1)<1/58 | abs(ptsec(5,:)-m+1)<1/58);
  subxpt = xpt(:,abs(ptsec(5,:)-m)<1/58 | abs(ptsec(5,:)-m-1)<1/58 | abs(ptsec(5,:)-m+1)<1/58);
  subdps = dps(abs(ptsec(5,:)-m)<1/58 | abs(ptsec(5,:)-m-1)<1/58 | abs(ptsec(5,:)-m+1)<1/58,:);
  C = cov(subY(1:4,:)'-subxpt');
  egs = sort(eigs(C));
  figure(4)
  hold on
  plot(m, std(subdps(:,1)),'ko','MarkerFaceColor', 'k')
  plot(m, mean(subdps(:,1)), 'bo', 'MarkerFaceColor', 'b')
  figure(5)
  hold on
  plot(m, std(subdps(:,2)),'ko','MarkerFaceColor', 'k')
  plot(m, mean(subdps(:,2)), 'bo', 'MarkerFaceColor', 'b')
end

figure(4)
set(gca, 'Fontsize', 14, 'LineWidth', 2)
xlabel('$\psi$','Interpreter','Latex','Fontsize',18)
ylabel('$\left(x-\gamma(\psi(x),0.5)\right)^\mathsf{T}e_1(\psi(x),0.5)$', ...
  'Interpreter', 'Latex', 'Fontsize', 18)

figure(5)
set(gca,'Fontsize',14, 'LineWidth', 2)
xlabel('$\psi$','Interpreter','Latex','Fontsize',18)
ylabel('$\left(x-\gamma(\psi(x),0.5)\right)^\mathsf{T}e_2(\psi(x),0.5)$', ...
  'Interpreter', 'Latex', 'Fontsize', 18)
