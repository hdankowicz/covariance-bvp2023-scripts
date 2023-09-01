%% Analysis of noise-perturbed quasiperiodic invariant torus
%
% Script generating Figure 6 in Z. Ahsan, H. Dankowicz, C. Kuehn,
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

%% Two glued instances of Tori with v(phi,tau)=u(om1*tau,phi+om2*tau)
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

prob = coco_prob;
prob = coco_set(prob, 'coll', 'MXCL', false, 'NTST', 20);

% Poincare conditions
data = VanderPol_init_data(funcs{1}, X0, p0, 1/rho, N, nsegs, xdim, pdim);
prob = ode_isol2bvp(prob, 'torus1', coll, {'eps1' 'beta1' 'del1'}, ...
  @VanderPol_bc, @VanderPol_bc_du, @VanderPol_bc_dudu, data, ...
  @Poincare_update); % first instance
prob = ode_isol2bvp(prob, 'torus2', coll, {'eps2' 'beta2' 'del2'}, ...
  @VanderPol_bc, @VanderPol_bc_du, @VanderPol_bc_dudu, data, ...
  @Poincare_update); % second instance

x0_idx = zeros(xdim*nsegs,1); x1_idx = zeros(xdim*nsegs,1);
for i=1:nsegs
  fname = sprintf('torus1.bvp.seg%d.coll',i);
  [uidx,fdata] = coco_get_func_data(prob,fname,'uidx','data');
  x0_idx(xdim*(i-1)+1:xdim*(i-1)+xdim) = uidx(fdata.coll_seg.maps.x0_idx);
  x1_idx(xdim*(i-1)+1:xdim*(i-1)+xdim) = uidx(fdata.coll_seg.maps.x1_idx);
end

% Rigid rotation
prob = coco_add_func(prob,'rot_bc1', @VanderPol_rot_bc, ...
  @VanderPol_rot_bc_du, @VanderPol_rot_bc_dudu, data, 'zero', ...
  'uidx', [x0_idx; x1_idx], 'u0', 1/rho);
uidx = coco_get_func_data(prob, 'rot_bc1', 'uidx');
prob = coco_add_pars(prob, 'rho1', uidx(end), 'rho1');

x0_idx = zeros(xdim*nsegs,1); x1_idx = zeros(xdim*nsegs,1);
for i=1:nsegs
  fname = sprintf('torus2.bvp.seg%d.coll',i);
  [uidx,fdata] = coco_get_func_data(prob,fname,'uidx','data');
  x0_idx(xdim*(i-1)+1:xdim*(i-1)+xdim) = uidx(fdata.coll_seg.maps.x0_idx);
  x1_idx(xdim*(i-1)+1:xdim*(i-1)+xdim) = uidx(fdata.coll_seg.maps.x1_idx);
end

% Rigid rotation
prob = coco_add_func(prob,'rot_bc2', @VanderPol_rot_bc, ...
  @VanderPol_rot_bc_du, @VanderPol_rot_bc_dudu, data, 'zero', ...
  'uidx', [x0_idx; x1_idx], 'u0', 1/rho);
uidx = coco_get_func_data(prob, 'rot_bc2', 'uidx');
prob = coco_add_pars(prob, 'rho2', uidx(end), 'rho2');

% gluing the problem parameter of the two instances
[uidx1, fdata1] = coco_get_func_data(prob, 'torus1.bvp.seg1.coll', ...
  'uidx', 'data');
[uidx2, fdata2] = coco_get_func_data(prob, 'torus2.bvp.seg1.coll', ...
  'uidx', 'data');
prob = coco_add_glue(prob, 'glue', ...
  uidx1(fdata1.coll_seg.maps.p_idx(2)), ...
  uidx2(fdata2.coll_seg.maps.p_idx(2)));

prob = coco_set(prob, 'cont', 'h_max', 5, 'NPR', 100);
prob = coco_set(prob, 'corr', 'MaxStep', 2);
prob = coco_add_event(prob, 'UZ', 'beta1', 0.5);
coco(prob, 'two_tori', [], 1, ...
  {'beta1' 'del1' 'beta2' 'del2'}, {[] [1.5 2]});

%% Final Covariance Problem
run = 'two_tori';
lab = coco_bd_labs('two_tori', 'UZ');
prob = coco_prob;
prob = coco_set(prob, 'coll', 'MXCL', false);

% Setting up the collocation problem
prob = ode_bvp2bvp(prob, 'torus1', run, lab);
[uidx, fdata] = coco_get_func_data(prob, 'torus1.bvp.seg1.coll', ...
  'uidx', 'data');
prob = coco_add_pars(prob, 'torus1_T', uidx(fdata.coll_seg.maps.T_idx), ...
  'torus1_T','active');

prob = ode_bvp2bvp(prob, 'torus2', run, lab);
[uidx, fdata] = coco_get_func_data(prob, 'torus2.bvp.seg1.coll', ...
  'uidx', 'data');
prob = coco_add_pars(prob, 'torus2_T', uidx(fdata.coll_seg.maps.T_idx), ...
  'torus2_T', 'active');

prob = adjt_isol2bvp(prob, 'torus1');
[axidx, opt] = coco_get_adjt_data(prob, 'torus1.bvp.seg1.coll', ...
  'axidx', 'data');
prob = coco_add_adjt(prob, 'torus1_T', 'd.torus1_T', 'aidx', ...
  axidx(opt.data.coll_opt.T_idx), 'l0', nsegs);

prob = adjt_isol2bvp(prob, 'torus2');
[axidx, opt] = coco_get_adjt_data(prob, 'torus2.bvp.seg1.coll', ...
  'axidx', 'data');
prob = coco_add_adjt(prob, 'torus2_T', 'd.torus2_T', 'aidx', ...
  axidx(opt.data.coll_opt.T_idx));

x0_idx = zeros(xdim*nsegs,1); x1_idx = zeros(xdim*nsegs,1);
for i=1:nsegs
  fname = sprintf('torus1.bvp.seg%d.coll',i);
  [uidx,fdata] = coco_get_func_data(prob,fname,'uidx','data');
  x0_idx(xdim*(i-1)+1:xdim*(i-1)+xdim) = uidx(fdata.coll_seg.maps.x0_idx);
  x1_idx(xdim*(i-1)+1:xdim*(i-1)+xdim) = uidx(fdata.coll_seg.maps.x1_idx);
end

% Rigid rotation
prob = coco_add_func(prob,'rot_bc1', @VanderPol_rot_bc, ...
  @VanderPol_rot_bc_du, @VanderPol_rot_bc_dudu, data, 'zero', ...
  'uidx', [x0_idx; x1_idx], 'u0', 1/rho);
uidx = coco_get_func_data(prob, 'rot_bc1', 'uidx');
prob = coco_add_pars(prob, 'rho1', uidx(end), 'rho1');

prob = coco_add_adjt(prob, 'rho1', 'd.rho1');

% adding adjoint of the rotation boundary condition
x0_adjt_idx = zeros(nsegs*xdim,1); x1_adjt_idx = zeros(nsegs*xdim,1);
for i=1:nsegs
  fname = sprintf('torus1.bvp.seg%d.coll',i);
  [axidx,fdata] = coco_get_adjt_data(prob,fname, 'axidx', 'data');
  x0_adjt_idx(xdim*(i-1)+1:xdim*(i-1)+xdim) = axidx(fdata.coll_opt.x0_idx);
  x1_adjt_idx(xdim*(i-1)+1:xdim*(i-1)+xdim) = axidx(fdata.coll_opt.x1_idx);
end

axidx_rho = coco_get_adjt_data(prob, 'rho1', 'axidx');
prob = coco_add_adjt(prob, 'rot_bc1', 'aidx', ...
  [x0_adjt_idx; x1_adjt_idx; axidx_rho]);

x0_idx = zeros(xdim*nsegs,1); x1_idx = zeros(xdim*nsegs,1);
for i=1:nsegs
  fname = sprintf('torus2.bvp.seg%d.coll',i);
  [uidx,fdata] = coco_get_func_data(prob,fname,'uidx','data');
  x0_idx(xdim*(i-1)+1:xdim*(i-1)+xdim) = uidx(fdata.coll_seg.maps.x0_idx);
  x1_idx(xdim*(i-1)+1:xdim*(i-1)+xdim) = uidx(fdata.coll_seg.maps.x1_idx);
end

% Rigid rotation
prob = coco_add_func(prob,'rot_bc2', @VanderPol_rot_bc, ...
  @VanderPol_rot_bc_du, @VanderPol_rot_bc_dudu, data, 'zero', ...
  'uidx', [x0_idx; x1_idx], 'u0', 1/rho);
uidx = coco_get_func_data(prob, 'rot_bc2', 'uidx');
prob = coco_add_pars(prob, 'rho2', uidx(end), 'rho2');

prob = coco_add_adjt(prob, 'rho2', 'd.rho2', 'l0', -nsegs);

% adding adjoint of the rotation boundary condition
x0_adjt_idx = zeros(nsegs*xdim,1); x1_adjt_idx = zeros(nsegs*xdim,1);
for i=1:nsegs
  fname = sprintf('torus2.bvp.seg%d.coll',i);
  [axidx,fdata] = coco_get_adjt_data(prob,fname, 'axidx', 'data');
  x0_adjt_idx(xdim*(i-1)+1:xdim*(i-1)+xdim) = axidx(fdata.coll_opt.x0_idx);
  x1_adjt_idx(xdim*(i-1)+1:xdim*(i-1)+xdim) = axidx(fdata.coll_opt.x1_idx);
end

axidx_rho = coco_get_adjt_data(prob, 'rho2', 'axidx');
prob = coco_add_adjt(prob, 'rot_bc2', 'aidx', ...
  [x0_adjt_idx; x1_adjt_idx; axidx_rho]);

% gluing the problem parameter of the two instances
[uidx1, fdata1] = coco_get_func_data(prob, 'torus1.bvp.seg1.coll', ...
  'uidx', 'data');
[uidx2, fdata2] = coco_get_func_data(prob, 'torus2.bvp.seg1.coll', ...
  'uidx', 'data');
prob = coco_add_glue(prob, 'glue', ...
  uidx1(fdata1.coll_seg.maps.p_idx(2)), ...
  uidx2(fdata2.coll_seg.maps.p_idx(2)));

data.fhan      = vanderpol('');
data.dfdxhan   = vanderpol('x');
data.Fnoisehan = @VanderPol_Fnoise;
data.uidx      = coco_get_func_data(prob, 'efunc', 'uidx');

prob = coco_add_slot(prob, '', @covariance, data, 'bddat');

prob = coco_set(prob, 'corr', 'MaxStep', 2);
coco(prob, 'two_tori_with_adjoints', [], 0, {'del1' 'beta2' 'del2' ...
  'd.eps1' 'd.beta1' 'd.del1' 'd.eps2' 'd.beta2' 'd.del2'});

%% Graphing Results

chart = coco_read_solution('rho1', 'two_tori_with_adjoints', 1, 'chart');
rho = chart.x;
data.rho  = rho;

[~, cdata] = coll_read_solution('torus1.bvp.seg1', ...
  'two_tori_with_adjoints', 1);
NTST = cdata.coll_seg.maps.NTST;
NCOL = cdata.coll_seg.int.NCOL;
data.NTST = NTST;
data.NCOL = NCOL;

covar = coco_bd_col('two_tori_with_adjoints', 'covariance');
Cbp_rep = zeros(NTST*(NCOL+1)*16,nsegs);
for i=1:nsegs
  Cbp_rep(:,i) = reshape(covar(:,:,:,i),NTST*(NCOL+1)*16,1);
end

data.ptmesh = Cbp_rep;
phi = 0:0.01:1;
tau = 0:0.01:1;
vals = zeros(101,101,2);
for i=1:101
  Crv = reshape(findCovars(data, [phi; tau(i)*ones(1,101)]'), [4,4,101]);
  for j=1:101
    eig = sort(eigs(Crv(:,:,j)));
    vals(i,j,:) = 0.1*sqrt(eig(3:4));
  end
end

figure
surf(phi, tau, vals(:,:,2))
grid on
box on
alpha 0.7
set(gca, 'Fontsize', 14, 'LineWidth', 2)
xlabel('$\phi$','Interpreter','Latex','Fontsize',20)
ylabel('$t$','Interpreter','Latex','Fontsize',20)
title('First largest eigenvalue','Interpreter','Latex','Fontsize',18)
view(-25,45)

figure
surf(phi, tau, vals(:,:,1))
grid on
box on
alpha 0.7
set(gca, 'Fontsize', 14, 'LineWidth', 2)
xlabel('$\phi$','Interpreter','Latex','Fontsize',20)
ylabel('$t$','Interpreter','Latex','Fontsize',20)
title('Second largest eigenvalue','Interpreter','Latex','Fontsize',18)
view(-25,45)

%% Covariance slot function

function [data, res] = covariance(prob, data, command, varargin)
%COVARIANCE   Append covariance data to BD.

res = {};
switch command
  
  case 'init'
    res   = { 'alpha' 'covariance' };
    
  case 'data'
    chart = varargin{1};
    fdata = coco_get_func_data(prob, 'torus1.bvp.seg1.coll', 'data');
    maps = fdata.coll_seg.maps;
    NTST = maps.NTST;
    NCOL = fdata.coll_seg.int.NCOL;
    xdim = fdata.xdim;
    N    = data.N;
    nsegs = data.nsegs;
    In = eye(xdim); In2 = eye(xdim^2);
    
    xbp   = cell(1,nsegs);
    xcn   = cell(1,nsegs);
    for i=1:nsegs
      uidx_i = coco_get_func_data(prob, ...
        sprintf('torus1.bvp.seg%d.coll', i), 'uidx');
      xbp{i} = chart.x(uidx_i(maps.xbp_idx));
      xcn{i} = reshape(maps.W*xbp{i}, maps.x_shp);
    end
    Xcn = maps.W*cell2mat(xbp);
    
    coeff  = cell(1,nsegs);
    alpha1 = [];
    alpha2a = [];
    alpha2b = [];
    alpha3 = [];
    rhs    = [];
    for i=1:nsegs
      [uidx_i] = coco_get_func_data(prob, ...
        sprintf('torus1.bvp.seg%d.coll', i), 'uidx');
      lidx_T_i = coco_get_adjt_data(prob, ...
        sprintf('torus1.bvp.seg%d.coll', i), 'afidx');
      lidx_rho_i = coco_get_adjt_data(prob, ...
        sprintf('torus2.bvp.seg%d.coll', i), 'afidx');
      
      T = chart.x(uidx_i(maps.T_idx));
      p = chart.x(uidx_i(maps.p_idx));
      lbp_T   = chart.x(data.uidx(end)+lidx_T_i);
      lbp_rho = chart.x(data.uidx(end)+lidx_rho_i);
      
      lcn_T   = reshape(maps.W*lbp_T, maps.x_shp);
      lcn_rho = reshape(maps.W*lbp_rho, maps.x_shp);
      
      fcn  = data.fhan(xcn{i},repmat(p,maps.p_rep));
      dfdx = data.dfdxhan(xcn{i},repmat(p,maps.p_rep));
      
      rows = reshape(1:NTST*NCOL*xdim^2, [xdim, NTST*NCOL*xdim]);
      rows = repmat(rows, xdim, 1);
      cols = repmat(1:NTST*NCOL*xdim^2, xdim, 1);
      vals = reshape(dfdx, [xdim^2, NTST*NCOL]);
      vals = repmat(vals, xdim, 1);
      InJ  = sparse(rows, cols, vals);
      JIn  = kron(sparse(maps.fdxrows, maps.fdxcols, dfdx(:)), In);
      
      W  = kron(maps.W, In);
      Wp = kron(maps.Wp, In);
      
      coeff{i} = [2*NTST/T*Wp-(InJ+JIn)*W; kron(maps.Q, In)];
      
      fv   = kron(fcn, In);
      rows = 1:NTST*NCOL*xdim;
      cols = repmat(1:NTST*NCOL, xdim, 1);
      fh   = sparse(rows, cols, fcn);
      fTf  = fv*fh;
      
      alpha1 = [alpha1; fTf(:); zeros((NTST-1)*xdim^2,1)]; %#ok<AGROW>
      
      phi = (i-1)/(2*N+1);  %===angle corresponding to the characteristic
      Th  = (1:N)*2*pi*phi;
      negkSinkCos = 2*pi*repmat(1:N,2,1).*[-sin(Th); cos(Th)];
      dvdphi = reshape(Xcn*data.F1d'*[0; negkSinkCos(:)], [xdim, NTST*NCOL]);
      
      dvdphiv = kron(dvdphi, In);
      dvdphih  = sparse(rows, cols, dvdphi);
      fTdvdphi = fv*dvdphih;
      dvdphiTf = dvdphiv*fh;
      dvdphiTdvdphi = dvdphiv*dvdphih;
      
      alpha2a = [alpha2a; fTdvdphi(:); zeros((NTST-1)*xdim^2,1)]; %#ok<AGROW>
      alpha2b = [alpha2b; dvdphiTf(:); zeros((NTST-1)*xdim^2,1)]; %#ok<AGROW>
      alpha3 = [alpha3; dvdphiTdvdphi(:); zeros((NTST-1)*xdim^2,1)]; %#ok<AGROW>
      
      %====Generating the projection matrix: Q
      rows = 1:NTST*NCOL*xdim;
      cols = repmat(1:NTST*NCOL, xdim, 1);
      DVDTAU = sparse(rows, cols, fcn);
      DVDPHI = sparse(rows, cols, dvdphi);
      LAMBDA_T   = sparse(cols, rows, lcn_T);
      LAMBDA_RHO = sparse(cols, rows, lcn_rho);
      
      Q = kron(eye(NTST*NCOL), In) - DVDTAU*LAMBDA_T - DVDPHI*LAMBDA_RHO;
      
      Fnoise = data.Fnoisehan(xcn{i}, repmat(p, maps.p_rep));
      Fnoise = sparse(maps.fdxrows, maps.fdxcols, Fnoise(:));
      Bcn    = Q*(Fnoise*Fnoise')*Q';
      
      rows = reshape(1:xdim*NTST*NCOL, [xdim, NTST*NCOL]);
      rows = repmat(rows, xdim, 1);
      cols = repmat(1:NTST*NCOL*xdim, xdim, 1);
      idx  = sub2ind(size(Bcn), rows(:), cols(:));
      
      rhs = [rhs; Bcn(idx); zeros((NTST-1)*xdim^2,1)]; %#ok<AGROW>
    end
    
    lidx_T_i = coco_get_adjt_data(prob, ...
      'torus1.bvp.seg1.coll', 'afidx');
    lidx_rho_i = coco_get_adjt_data(prob, ...
      'torus2.bvp.seg1.coll', 'afidx');
    
    lbp_T   = chart.x(data.uidx(end)+lidx_T_i);
    lbp_rho = chart.x(data.uidx(end)+lidx_rho_i);
    
    temp = [blkdiag(coeff{:}) -alpha1 -alpha2a -alpha2b -alpha3;
      kron(lbp_T(1:xdim)',lbp_T(1:xdim)'), ...
      zeros(1,(nsegs*NTST*(NCOL+1)-1)*xdim^2+4);
      kron(lbp_T(1:xdim)',lbp_rho(1:xdim)'), ...
      zeros(1,(nsegs*NTST*(NCOL+1)-1)*xdim^2+4);
      kron(lbp_rho(1:xdim)',lbp_T(1:xdim)'), ...
      zeros(1,(nsegs*NTST*(NCOL+1)-1)*xdim^2+4);
      kron(lbp_rho(1:xdim)',lbp_rho(1:xdim)'), ...
      zeros(1,(nsegs*NTST*(NCOL+1)-1)*xdim^2+4)];
    
    F_bc  = kron(data.F1d, In2);
    RF_bc = kron(data.R*data.F1d, In2);
    
    bc = zeros(nsegs*xdim^2,nsegs*NTST*(NCOL+1)*xdim^2);
    for i=1:nsegs
      bc(:,NTST*(NCOL+1)*xdim^2*(i-1)+(1:xdim^2)) = ...
        RF_bc(:,xdim^2*(i-1)+(1:xdim^2));
      bc(:,NTST*(NCOL+1)*xdim^2*i+(1-xdim^2:0)) = ...
        -F_bc(:,xdim^2*(i-1)+(1:xdim^2));
    end
    
    dump = [temp; bc, zeros(nsegs*xdim^2,4)] \ ...
      [rhs; zeros(nsegs*xdim^2+4,1)];
    
    res = { dump(end-3:end), ...
      reshape(dump(1:end-4), xdim, xdim, NTST*(NCOL+1), nsegs) };
end

end

% the covariance problem:
% cdot = J*C+C*J^T+Q*F*F^T*Q^T+alpha*f*f^T, c(0)=c(T), l(0)^T*c(0)*l(0)=0
% vectorization:
% cvecdot - (I_n ox J + J ox I_n)*cvec - alpha*(f ox f) = (Q*F) ox (Q*F),
% RF_bc*C0-F_bc*C1, (l_T(0)^T ox l_T(0)^T)*cvec(0)=0,
% (l_T(0)^T ox l_rho(0)^T)*cvec(0)=0, (l_rho(0)^T ox l_rho(0)^T)*cvec(0)=0

%% noise model
function F = VanderPol_Fnoise(x,p) %#ok<INUSD>

x1 = x(1,:);

F = zeros(4,4,numel(x1));

F(2,1,:) = 1;
F(4,2,:) = 1;

end

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
