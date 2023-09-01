%% Analysis of noise-perturbed quasiperiodic invariant torus
%
% Script generating Figure 4 in Z. Ahsan, H. Dankowicz, C. Kuehn,
% "Adjoint-Based Projections for Uncertainty Quantification near
% Stochastically Perturbed Limit Cycles and Tori"

%% Simulation
T  = 10000;
dt = 10^(-4);
t  = 0:dt:T;
N  = ceil(T/dt);

x1  = sqrt(2);
x2  = 0;
Xem = [x1 x2; zeros(N,2)];

% Generating the noise realization
% Storing it on disk to save memory
dW = sqrt(dt)*randn(N,1);
save('Noise_realization_quasiperiodic_ex1_10000_01.mat','dW')
% Simulation setup

sigma = 0.1;
omega = 1;
rho   = pi;
Omega = rho*omega;
for j=1:numel(dW)
  Fnoise = [x1*x2*dW(j); x2^2*dW(j)];
  r = 1+sqrt(x1^2+x2^2)*(cos(2*pi*t(j))-1);
  f = [-Omega*x2+x1*r; Omega*x1+x2*r];
  x1 = x1+dt*2*pi/omega*f(1)+sigma*sqrt(2*pi/omega)*Fnoise(1);
  x2 = x2+dt*2*pi/omega*f(2)+sigma*sqrt(2*pi/omega)*Fnoise(2);
  Xem(j+1,:) = [x1 x2];
end

save("stochastic_quasiperiodic_ex1_10000_01.mat","t","Xem")

%% Visualization

load("stochastic_quasiperiodic_ex1_10000_01.mat","t","Xem")

secn = 0;
Y=[];
c=0;
for i=1:length(t)
  tmod = t(i)-floor(t(i));
  if abs(tmod-secn)<1e-8
    c = c+1;
    Y(c,:) = Xem(i,:);
  end
end

figure
grid on
box on
hold on
plot(Y(:,1),Y(:,2),'r.','Markersize',4)

phi = 0:0.01:1;
r   = (1+omega^2)/(1+omega^2-cos(2*pi*secn)-omega*sin(2*pi*secn));
x   = r*cos(2*pi*(phi+rho*secn));
y   = r*sin(2*pi*(phi+rho*secn));

plot(x, y, 'Linewidth', 2)

lambda = sigma^2*(omega^2+1)^4*(1+Omega^2-cos(4*pi*phi) ...
  -Omega*sin(4*pi*phi))/4/omega^8/(1+Omega^2);

mag = sqrt(x.^2+y.^2);
for radius = -2:2
  if radius~=0
    plot(x+radius*sqrt(lambda).*(x./mag), ...
      y+radius*sqrt(lambda).*(y./mag), 'k', 'LineWidth', 2)
  else
    plot(x-radius*sqrt(lambda).*(x./mag), ...
      y-radius*sqrt(lambda).*(y./mag), 'k--', 'LineWidth', 2)
  end
end
axis equal
axis tight
hold off
set(gca, 'Linewidth', 2, 'Fontsize', 14)
xlabel('$x_1$', 'Interpreter', 'Latex', 'Fontsize', 20)
ylabel('$x_2$', 'Interpreter', 'Latex', 'Fontsize', 20)

figure
hold on
grid on
box on
angles = atan2(Y(:,2),Y(:,1));
angles(angles<0) = angles(angles<0)+2*pi;
plot(angles/(2*pi),sqrt(Y(:,1).^2+Y(:,2).^2)-r,'r.')
for radius=-2:2
  if radius~=0
    plot(phi, radius*sqrt(lambda), 'k', 'LineWidth', 2)
  else
    plot(phi, radius*sqrt(lambda), 'k--', 'LineWidth', 2)
  end
end
for m=0:5*pi/100:2*pi
  subY = Y(abs(angles-m)<5*pi/200 | abs(angles-m-2*pi)<5*pi/200 ...
    | abs(angles-m+2*pi)<5*pi/200,:);
  plot(m/(2*pi), std(sqrt(subY(:,1).^2+subY(:,2).^2)), ...
    'ko', 'MarkerFaceColor', 'k')
  plot(m/(2*pi), mean(sqrt(subY(:,1).^2+subY(:,2).^2))-2, ...
    'bo', 'MarkerFaceColor', 'b')
end
hold off
axis tight
set(gca, 'Linewidth', 2, 'Fontsize', 14)
xlabel('$\psi$', 'Interpreter', 'Latex', 'Fontsize', 20)
ylabel('$\left(x(t)-\gamma(\psi(t),\tau(t))\right)^\mathsf{T}e(\psi(t),\tau(t))$', ...
  'Interpreter', 'Latex', 'Fontsize', 20)
title('$t\in \{0,1,2,\ldots\}$', 'Interpreter', 'Latex', 'Fontsize', 18)
