%% Analysis of noise-perturbed limit cycle
%
% Script generating Figure 1 in Z. Ahsan, H. Dankowicz, C. Kuehn,
% "Adjoint-Based Projections for Uncertainty Quantification near
% Stochastically Perturbed Limit Cycles and Tori"

%% Creating functions

if exist('sys_hopf', 'file')~=2
  syms x1 x2 mu
  f = [mu*x1-x2-x1*(x1^2+x2^2); x1+mu*x2-x2*(x1^2+x2^2)];
  hopf = sco_sym2funcs(f, {[x1; x2], mu}, {'x', 'p'}, 'filename', 'sys_hopf');
else
  hopf = sco_gen(@sys_hopf);
end

funcs = {hopf(''), hopf('x'), hopf('p'), hopf({'x','x'}), ...
  hopf({'x','p'}), hopf({'p','p'})};

%% Solving the periodic orbit along with the adjoint problem
t0 = (0:0.01:2*pi)';
x0 = [cos(t0) sin(t0)];

prob = coco_prob;
prob = coco_set(prob, 'coll', 'NTST', 10);
prob = ode_isol2po(prob, '', funcs{:}, t0, x0, 'mu', 1);

prob = adjt_isol2po(prob, '');
Tidx = coco_get_adjt_data(prob, 'po.period', 'afidx');
fcn  = @(f) @(p,d,u,l,v) deal(d, f(u,l,v));
prob = coco_add_comp(prob, 'fun1', fcn(@(u,l,v) l-1), [], 'zero', ...
  'lidx', Tidx);

prob = coco_set(prob, 'corr', 'MaxStep', 2);
coco(prob, 'run1', [], 0, {'d.po.period','d.mu'});

%% Solving the periodic orbit and covariance problem
prob = coco_prob;
prob = ode_po2po(prob,'','run1', 1);

prob = adjt_po2po(prob,'','run1', 1);
Tidx = coco_get_adjt_data(prob, 'po.period', 'afidx');
fcn  = @(f) @(p,d,u,l,v) deal(d, f(u,l,v));
prob = coco_add_comp(prob, 'fun1', fcn(@(u,l,v) l-1), [], 'zero', ...
  'lidx', Tidx);

data = struct();
data.fhan      = hopf('');
data.dfdxhan   = hopf('x');
data.Fnoisehan = @Hopf_noise;
data.fbid      = 'po.orb.coll';
data.uidx      = coco_get_func_data(prob, 'efunc', 'uidx');

prob = coco_add_slot(prob, '', @covariance, data, 'bddat');

coco(prob, 'run2', [], 0, {'d.po.period','d.mu'});

%% Comparing computed eigenvalue to theoretical prediction
covar = coco_bd_col('run2','covariance');
las = zeros(2,50);
for i=1:50
  las(:,i) = eig(covar(:,:,i));
end
[~, data] = coll_read_solution('po.orb','run2',1);
tbp = data.coll_seg.mesh.tbp;

figure
hold on
grid on
box on
plot(tbp, 1/40*(5-4*cos(4*pi*tbp)-2*sin(4*pi*tbp)), ...
  'b', 'LineWidth', 2, 'DisplayName', 'Analytical')
plot(tbp, max(las), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'COCO')
set(gca, 'FontSize', 14, 'Linewidth', 2);
xlabel('$\tau$','Interpreter','Latex','Fontsize',20)
legend('Position', [0.7 0.2 0.1 0.1], 'FontSize', 15);
legend('boxoff')
axis([0 1 0 0.25])
set(gcf,'position',[0,200,430,310])

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
function F = Hopf_noise(x,p) %#ok<INUSD>

x1 = x(1,:); 
x2 = x(2,:);

F = zeros(2,2,numel(x1));

F(1,1,:) = x1.*x2;
F(2,1,:) = x2.^2;

end
