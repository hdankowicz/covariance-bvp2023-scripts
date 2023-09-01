%% Analysis of noise-perturbed quasiperiodic invariant torus
%
% Script generating Figure 5 in Z. Ahsan, H. Dankowicz, C. Kuehn,
% "Adjoint-Based Projections for Uncertainty Quantification near
% Stochastically Perturbed Limit Cycles and Tori"

%% Creating functions

if exist('sys_vanderpol', 'file')~=2
  syms x1 x2 x3 x4 eps beta del
  f = [x2; -eps*(x1^2-1)*x2-x1+beta*(x3-x1); ...
    x4; -eps*(x3^2-1)*x4-(1+del)*x3+beta*(x1-x3)];
  vanderpol = sco_sym2funcs(f, {[x1; x2; x3; x4], [eps; beta; del]}, ...
    {'x', 'p'}, 'filename', 'sys_vanderpol');
else
  vanderpol = sco_gen(@sys_vanderpol);
end

funcs = {vanderpol(''), vanderpol('x'), vanderpol('p'), ...
  vanderpol({'x','x'}), vanderpol({'x','p'}), vanderpol({'p','p'})};

%% Single torus with v(phi,tau)=u(om1*tau,phi+om2*tau)
% in this case, the rotation number is om2/om1 = 1/rho

om1 = 1;
n   = 62;
rho = n*sqrt(2)/140;
om2 = om1/rho;

eps = 0.5; beta = 0; del = om2^2-1;
p0  = [eps; beta; del];

N = 14; nsegs = 2*N+1; xdim = 4; pdim = 3;
tau = 2*pi/om1*linspace(0,1,10*nsegs);
th1 = om1*tau(:);

coll = cell(1,nsegs);
X0   = zeros(xdim,nsegs);
for i=1:nsegs
  phi = 2*pi*(i-1)/nsegs;
  th2 = phi + om2*tau(:);
  up  = [2*sin(th1(:)) 2*om1*cos(th1(:)) 2*sin(th2(:)) 2*om2*cos(th2(:))];
  coll{i} = [funcs, {tau}, {up}, {[eps beta del]}];
  X0(:,i) = up(1,:)';
end

prob = coco_prob();
prob = coco_set(prob, 'coll', 'MXCL', false, 'NTST', 20);

% Poincare conditions
data = VanderPol_init_data(funcs{1}, X0, p0, 1/rho, N, nsegs, xdim, pdim);
prob = ode_isol2bvp(prob, 'torus', coll, {'eps' 'beta' 'del'}, ...
  @VanderPol_bc, @VanderPol_bc_du, @VanderPol_bc_dudu, data, ...
  @Poincare_update);

x0_idx = zeros(xdim*nsegs,1); x1_idx = zeros(xdim*nsegs,1);
for i=1:nsegs
  fname = sprintf('torus.bvp.seg%d.coll',i);
  [uidx,fdata] = coco_get_func_data(prob,fname, 'uidx', 'data');
  x0_idx(xdim*(i-1)+1:xdim*(i-1)+xdim) = uidx(fdata.coll_seg.maps.x0_idx);
  x1_idx(xdim*(i-1)+1:xdim*(i-1)+xdim) = uidx(fdata.coll_seg.maps.x1_idx);
end

% Rigid rotation
prob = coco_add_func(prob,'rot_bc', @VanderPol_rot_bc, ...
  @VanderPol_rot_bc_du, @VanderPol_rot_bc_dudu, data, 'zero', ...
  'uidx', [x0_idx; x1_idx], 'u0', 1/rho);
uidx = coco_get_func_data(prob, 'rot_bc', 'uidx');
prob = coco_add_pars(prob, 'rho', uidx(end), 'rho');

prob = coco_set(prob, 'cont', 'h_max', 5, 'NPR', 100);
prob = coco_set(prob, 'corr', 'MaxStep', 2);
prob = coco_add_event(prob, 'UZ', 'beta', 0.5);
coco(prob, 'torus', [], 1, {'beta' 'del'}, {[] [1.5 2]});

%% Visualization
figure
thm = struct();
thm.lspec = {{'b', 'LineWidth', 2}};
coco_plot_bd(thm, 'torus', 'del', 'beta')
grid on
box on
set(gca, 'LineWidth', 2, 'Fontsize', 14)
xlabel('$\delta$', 'Interpreter', 'Latex', 'Fontsize', 20)
ylabel('$\beta$',  'Interpreter', 'Latex', 'Fontsize', 20)

figure
coco_plot_sol('torus', 4, 'torus', 'x', 'x', 'x')
box on
grid on
axis tight
view(-100,15)
set(gca, 'LineWidth', 2, 'Fontsize',14)
xlabel('$x_1$','Interpreter','Latex','Fontsize',20)
ylabel('$x_2$','Interpreter','Latex','Fontsize',20)
zlabel('$x_3$','Interpreter','Latex','Fontsize',20)

%% boundary conditions
function fbc = VanderPol_bc(data, T, x0, ~, ~)

xdim = data.xdim;

yphase1 = data.dvdtau0'*(x0(1:xdim,1)-data.x0(:,1));
yphase2 = data.dvdphi0'*(x0(1:xdim,1)-data.x0(:,1));

fbc = [T(1)-T(2:end); yphase1; yphase2];

end

function Jbc = VanderPol_bc_du(data, ~, ~, ~, ~)

nsegs = data.nsegs;
xdim  = data.xdim;
pdim  = data.pdim;

Jbc = zeros(nsegs+1, nsegs+2*nsegs*xdim+pdim);

for i=1:nsegs-1
  Jbc(i,1)   =  1;
  Jbc(i,i+1) = -1;
end

Jbc(nsegs,   nsegs+1:nsegs+xdim) = data.dvdtau0';
Jbc(nsegs+1, nsegs+1:nsegs+xdim) = data.dvdphi0';

end

function Jbc = VanderPol_bc_dudu(data, ~, ~, ~, ~, ~)

nsegs = data.nsegs;
xdim  = data.xdim;
pdim  = data.pdim;

Jbc = zeros(nsegs+1, nsegs+2*nsegs*xdim+pdim, nsegs+2*nsegs*xdim+pdim);

end

function [data,fbc] = VanderPol_rot_bc(~,data,u)

nsegs = data.nsegs;
xdim  = data.xdim;

x0 = u(1:nsegs*xdim,1);
x1 = u(nsegs*xdim+1:2*nsegs*xdim,1);

fbc = data.RF*x0-data.F*x1;

end

function [data, Jrot] = VanderPol_rot_bc_du(~,data,u)

nsegs = data.nsegs;
xdim  = data.xdim;

x0   = u(1:nsegs*xdim,1);
Jrot = zeros(nsegs*xdim,2*nsegs*xdim+1);
Jrot(1:nsegs*xdim,1:nsegs*xdim) = data.RF;
Jrot(1:nsegs*xdim,nsegs*xdim+1:2*nsegs*xdim) = -data.F;
Jrot(1:nsegs*xdim,end) = data.dRdrhoF*x0;

end

function [data, Jbc] = VanderPol_rot_bc_dudu(~, data, u)

nsegs = data.nsegs;
xdim  = data.xdim;

x0  = u(1:nsegs*xdim,1);
Jbc = zeros(nsegs*xdim,2*nsegs*xdim+1,2*nsegs*xdim+1);
Jbc(1:nsegs*xdim,1:nsegs*xdim,end) = data.dRdrhoF;
for j=1:nsegs*xdim
  dx0 = zeros(nsegs*xdim,1);
  dx0(j) = 1;
  Jbc(1:nsegs*xdim,end,j) = data.dRdrhoF*dx0;
end
Jbc(1:nsegs*xdim,end,end) = data.dRdrhodrhoF*x0;

end

%% initialization of data
function data = VanderPol_init_data(fhan, X0, p0, rho, N, nsegs, xdim, pdim)

data       = struct();
data.N     = N;
data.nsegs = nsegs;
data.xdim  = xdim;
data.pdim  = pdim;
data.fhan  = fhan;

Th  = 2*pi*(0:2*N)/(2*N+1);
Th  = kron(1:N, Th');
data.F1d = [ones(2*N+1,1) 2*reshape([cos(Th);sin(Th)], ...
  [2*N+1 2*N])]'/(2*N+1);
data.cfs = [zeros(1,N); 2*pi*(1:N)];

data = Poincare_update(data, [], X0(:), [], p0);

Th  = (1:N)*2*pi*rho;
SIN = [ zeros(size(Th)) ; sin(Th) ];
R   = diag([1 kron(cos(Th), [1 1])]);
data.R   = R  + diag(SIN(:), +1)- diag(SIN(:), -1);

data.dRdrho     = zeros(2*N+1,2*N+1);
data.dRdrhodrho = zeros(2*N+1,2*N+1);
for j=1:N
  Mat1 = [-2*pi*j*sin(2*pi*j*rho) 2*pi*j*cos(2*pi*j*rho);...
    -2*pi*j*cos(2*pi*j*rho) -2*pi*j*sin(2*pi*j*rho)];
  Mat2 = [-(2*pi*j)^2*cos(2*pi*j*rho) -(2*pi*j)^2*sin(2*pi*j*rho);...
    (2*pi*j)^2*sin(2*pi*j*rho) -(2*pi*j)^2*cos(2*pi*j*rho)];
  data.dRdrho(2*j:2*j+1,2*j:2*j+1) = Mat1;
  data.dRdrhodrho(2*j:2*j+1,2*j:2*j+1) = Mat2;
end

data.RF          = kron(data.R*data.F1d, eye(xdim));
data.dRdrhoF     = kron(data.dRdrho*data.F1d, eye(xdim));
data.dRdrhodrhoF = kron(data.dRdrhodrho*data.F1d, eye(xdim));
data.F           = kron(data.F1d, eye(xdim));

end

%% Poincare section update
function data = Poincare_update(data,~,x0,~,p)

data.x0 = reshape(x0, data.xdim, data.nsegs);
data.p  = p;

data.dvdtau0 = data.fhan(data.x0(:,1), data.p);
data.dvdphi0 = data.x0*data.F1d'*[0; data.cfs(:)];

end
