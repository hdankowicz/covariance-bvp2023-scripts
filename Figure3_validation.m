%% Analysis of noise-perturbed limit cycle
%
% Script validating encoding of covariance boundary-value problem
% corresponding to Figure 3 in Z. Ahsan, H. Dankowicz, C. Kuehn,
% "Adjoint-Based Projections for Uncertainty Quantification near
% Stochastically Perturbed Limit Cycles and Tori"

%% Creating functions

if exist('sys_linosc', 'file')~=2
  syms x1 x2 x3 p
  f = [x2; -p*x2-x1+2*cos(x3); 1];
  linosc = sco_sym2funcs(f, {[x1; x2; x3] p}, {'x' 'p'}, 'filename', 'sys_linosc');
else
  linosc = sco_gen(@sys_linosc);
end

funcs = {linosc(''), linosc('x'), linosc('p'), linosc({'x','x'}), ...
  linosc({'x','p'}), linosc({'p','p'})};

%% Solving the periodic orbit along with the adjoint problem
t0 = (0:0.01:2*pi)';
x0 = [sin(t0) cos(t0) t0];

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 10);
prob = ode_isol2bvp(prob, '', funcs{:}, t0, x0, 2, 'p', ...
  @linosc_bc, @linosc_bc_du, @linosc_bc_dudu);

[data, uidx] = coco_get_func_data(prob, 'bvp.seg1.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_pars(prob, 'period', uidx(maps.T_idx), 'period');

prob = adjt_isol2bvp(prob, '');
[data, axidx] = coco_get_adjt_data(prob, 'bvp.seg1.coll', 'data', 'axidx');
opt = data.coll_opt;
prob = coco_add_adjt(prob, 'period', 'd.period', 'aidx', axidx(opt.T_idx));
prob = coco_set_parival(prob, 'd.period', 1);

coco(prob, 'run1', [], 0, {'period', 'd.p'});

%% Solving the periodic orbit and covariance problem
prob = coco_prob();
prob = ode_bvp2bvp(prob,'','run1', 1);

[data, uidx] = coco_get_func_data(prob, 'bvp.seg1.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_pars(prob, 'period', uidx(maps.T_idx), 'period');

prob = adjt_bvp2bvp(prob,'','run1', 1);
[data, axidx] = coco_get_adjt_data(prob, 'bvp.seg1.coll', 'data', 'axidx');
opt = data.coll_opt;
prob = coco_add_adjt(prob, 'period', 'd.period', 'aidx', axidx(opt.T_idx));
prob = coco_set_parival(prob, 'd.period', 1);

data = struct();
data.fhan      = linosc('');
data.dfdxhan   = linosc('x');
data.Fnoisehan = @linosc_noise;
data.fbid      = 'bvp.seg1.coll';
data.uidx      = coco_get_func_data(prob, 'efunc', 'uidx');

prob = coco_add_slot(prob, '', @covariance, data, 'bddat');

coco(prob, 'run2', [], 0, {'period', 'd.p'});

%% Comparing computed eigenvalue to theoretical prediction
covar = coco_bd_col('run2','covariance');
las = zeros(3,50);
for i=1:50
  las(:,i) = eig(covar(:,:,i));
end
[~, data] = coll_read_solution('bvp.seg1','run2',1);
tbp = data.coll_seg.mesh.tbp;

figure
hold on
grid on
box on
set(gca, 'FontSize', 14);
plot(tbp, 1/(32)*(4-cos(4*pi*tbp)-sin(4*pi*tbp) ...
  -sqrt(3+2*cos(8*pi*tbp)+sin(8*pi*tbp))), ...
  'g', 'LineWidth', 2, 'DisplayName', 'Analytical')
plot(tbp, 1/(32)*(4-cos(4*pi*tbp)-sin(4*pi*tbp) ...
  +sqrt(3+2*cos(8*pi*tbp)+sin(8*pi*tbp))), ...
  'g', 'LineWidth', 2, 'DisplayName', 'Analytical')
plot(tbp, maxk(las,2), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'COCO')
xlabel('$\tau$','Interpreter','Latex','Fontsize',18)
set(gcf,'position',[0,200,430,310])
set(gca,'Linewidth',2)
hold off

%% Covariance slot function

function [data, res] = covariance(prob, data, command, varargin)
%COVARIANCE   Append covariance data to BD.

res = {};
switch command
  
  case 'init'
    res = { 'alpha', 'covariance' };
    
  case 'data'
    chart = varargin{1};
    lidx  = coco_get_adjt_data(prob, data.fbid, 'afidx');
    [uidx, fdata] = coco_get_func_data(prob, data.fbid, 'uidx', 'data');
    
    maps = fdata.coll_seg.maps;    
    NTST = maps.NTST;
    NCOL = fdata.coll_seg.int.NCOL;
    xdim = fdata.xdim;
    
    In = eye(xdim); In2 = eye(xdim^2);
    
    xbp = chart.x(uidx(maps.xbp_idx));
    T   = chart.x(uidx(maps.T_idx));
    p   = chart.x(uidx(maps.p_idx));
    lbp = chart.x(data.uidx(end)+lidx);
    
    xcn = reshape(maps.W*xbp, maps.x_shp);
    lcn = reshape(maps.W*lbp, maps.x_shp);

    dfdx = data.dfdxhan(xcn,repmat(p,maps.p_rep));
    rows = reshape(1:NTST*NCOL*xdim^2, [xdim, NTST*NCOL*xdim]);
    rows = repmat(rows, xdim, 1);
    cols = repmat(1:NTST*NCOL*xdim^2, xdim, 1);
    vals = reshape(dfdx, [xdim^2, NTST*NCOL]);
    vals = repmat(vals, xdim, 1);
    InJ  = sparse(rows, cols, vals);
    JIn  = kron(sparse(maps.fdxrows, maps.fdxcols, dfdx(:)), In);
    
    W  = kron(maps.W, In);
    Wp = kron(maps.Wp, In);
        
    fcn  = data.fhan(xcn,repmat(p,maps.p_rep));
    rows = 1:NTST*NCOL*xdim;
    cols = repmat(1:NTST*NCOL, xdim, 1);
    fh   = sparse(rows, cols, fcn);
    ld   = sparse(cols, rows, lcn);
    
    Q = kron(eye(NTST*NCOL), In) - fh*ld;
    
    fTf = kron(fcn, In)*fh;
    
    temp = [2*NTST/T*Wp-(InJ+JIn)*W -fTf(:);
      kron(maps.Q, In) zeros((NTST-1)*xdim^2,1);
      In2 zeros(xdim^2,(NTST*(NCOL+1)-2)*xdim^2) -In2 zeros(xdim^2,1);
      kron(lbp(1:xdim)',lbp(1:xdim)') zeros(1,(NTST*(NCOL+1)-1)*xdim^2+1)];
    
    Fnoise = data.Fnoisehan(xcn, repmat(p, maps.p_rep));
    Fnoise = sparse(maps.fdxrows, maps.fdxcols, Fnoise(:));
    Bcn    = Q*(Fnoise*Fnoise')*Q';
    
    rows = reshape(1:xdim*NTST*NCOL, [xdim, NTST*NCOL]);
    rows = repmat(rows, xdim, 1);
    cols = repmat(1:NTST*NCOL*xdim, xdim, 1);
    idx  = sub2ind(size(Bcn), rows(:), cols(:));
    
    dump = temp \ [Bcn(idx); zeros(NTST*xdim^2+1,1)];
    
    res = { dump(end), reshape(dump(1:end-1), xdim, xdim, NTST*(NCOL+1)) };
    
end

end

% the covariance problem:
% cdot = J*C+C*J^T+Q*F*F^T*Q^T+alpha*f*f^T, c(0)=c(T), l(0)^T*c(0)*l(0)=0
% vectorization:
% cvecdot - (I_n ox J + J ox I_n)*cvec - alpha*(f ox f) = (Q*F) ox (Q*F),
% cvec(0)=cvec(T), (l(0)^T ox l(0)^T)*cvec(0)=0

%% noise model
function F = linosc_noise(x,p) %#ok<INUSD>

x1 = x(1,:);

F = zeros(3,3,numel(x1));

F(2,1,:) = x1;

end

%% boundary conditions
function y = linosc_bc(~, ~, x0, x1, ~)
y = [x1(1:2)-x0(1:2); x1(3)-x0(3)-2*pi; x0(3)];
end

function J = linosc_bc_du(~, ~, ~, ~, ~)
J = zeros(4,8);
J(1,2) = -1;
J(1,5) = 1;
J(2,3) = -1;
J(2,6) = 1;
J(3,4) = -1;
J(3,7) = 1;
J(4,4) = 1;
end

function dJ = linosc_bc_dudu(~, ~, ~, ~, ~)
dJ = zeros(4,8,8);
end
