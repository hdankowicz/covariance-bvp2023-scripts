%% Analysis of noise-perturbed limit cycle
%
% Script generating Figure 2 in Z. Ahsan, H. Dankowicz, C. Kuehn,
% "Adjoint-Based Projections for Uncertainty Quantification near
% Stochastically Perturbed Limit Cycles and Tori"

%% Simulation
T  = 500;
dt = 10^(-4);
t  = 0:dt:T;
N  = ceil(T/dt);

x1  = 0;   
x2  = 1;
Xem = [x1 x2; zeros(N,2)];

% Generating the noise realization
% Storing it on disk to save memory
dW = sqrt(dt)*randn(N,1);
save('Noise_realization_periodic_ex1_500_01.mat','dW')
% Simulation setup

sigma = 0.1;
for j=1:numel(dW)
  Fnoise = [x1*x2*dW(j); x2^2*dW(j)];
  x1 = x1+dt*(x1-x2-x1*(x1^2+x2^2))*2*pi+sigma*sqrt(2*pi)*Fnoise(1);
  x2 = x2+dt*(x1+x2-x2*(x1^2+x2^2))*2*pi+sigma*sqrt(2*pi)*Fnoise(2);
  Xem(j+1,:) = [x1 x2];
end

save("stochastic_new_periodic_ex1_500_01.mat","t","Xem")

%% Visualization

load("stochastic_new_periodic_ex1_500_01.mat","t","Xem")

phi = 0:0.01:1;
lambda = sigma^2*(5-4*cos(4*pi*phi)-2*sin(4*pi*phi))/40;

figure
grid on
box on
hold on
plot(Xem(:,1), Xem(:,2), 'r.', 'Markersize', 1)
for radius=-2:2
  if radius~=0
    plot((1+radius*sqrt(lambda)).*cos(2*pi*phi), ...
      (1+radius*sqrt(lambda)).*sin(2*pi*phi), 'k', 'LineWidth', 2);
  else
    plot((1+radius*sqrt(lambda)).*cos(2*pi*phi), ...
      (1+radius*sqrt(lambda)).*sin(2*pi*phi), 'k--', 'LineWidth', 2);
  end
end
hold off
axis equal
axis([-1.3 1.3 -1.3 1.3])
set(gca, 'Linewidth', 2, 'Fontsize', 14)
xlabel('$x_1$', 'Interpreter', 'Latex', 'Fontsize', 20)
ylabel('$x_2$', 'Interpreter', 'Latex', 'Fontsize', 20)

figure
grid on
box on
hold on
angles = atan2(Xem(1:1:end,2),Xem(1:1:end,1));
angles(angles<0) = angles(angles<0)+2*pi;
plot(angles(1:1:end)/(2*pi), ...
  sqrt(Xem(1:1:end,1).^2+Xem(1:1:end,2).^2)-1, 'r.', 'MarkerSize', 0.5)
for radius=-2:2
  if radius~=0
    plot(phi,radius*sqrt(lambda), 'k', 'LineWidth', 2)
  else
    plot(phi,radius*sqrt(lambda), 'k', 'LineWidth', 2)
  end
end
for m=0:5*pi/100:2*pi
  subY = Xem(abs(angles-m)<5*pi/200 | abs(angles-m-2*pi)<5*pi/200 ...
    | abs(angles-m+2*pi)<5*pi/200,:);
  plot(m/(2*pi), std(sqrt(subY(:,1).^2+subY(:,2).^2)), 'ko', ...
    'MarkerFaceColor','k')
  plot(m/(2*pi), mean(sqrt(subY(:,1).^2+subY(:,2).^2))-1, 'bo', ...
    'MarkerFaceColor','b')
end
hold off
axis([0 1 -0.2 0.25])
set(gca, 'Linewidth', 2, 'Fontsize', 14)
xlabel('$\tau$', 'Interpreter', 'Latex', 'Fontsize', 20)
ylabel('$\left(x(t)-\gamma(\tau(t))\right)^\mathsf{T}e(\tau(t))$', ...
  'Interpreter', 'Latex', 'Fontsize', 20)
