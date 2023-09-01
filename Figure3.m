%% Analysis of noise-perturbed limit cycle
%
% Script generating Figure 3 in Z. Ahsan, H. Dankowicz, C. Kuehn,
% "Adjoint-Based Projections for Uncertainty Quantification near
% Stochastically Perturbed Limit Cycles and Tori"

%% Simulation
T  = 1000;
dt = 10^(-4);
t  = 0:dt:T;
N  = ceil(T/dt);

x1  = 0;
x2  = 1;
Xem = [x1 x2; zeros(N,2)];

% Generating the noise realization
% Storing it on disk to save memory
dW = sqrt(dt)*randn(N,1);
save('Noise_realization_periodic_ex2_1000_01.mat','dW')
% Simulation setup

sigma = 0.1;
for j=1:numel(dW)
  Fnoise = [0; x1*dW(j)];
  x1 = x1+dt*x2*2*pi+sigma*sqrt(2*pi)*Fnoise(1);
  x2 = x2+dt*(-2*x2-x1+2*cos(2*pi*t(j)))*2*pi+sigma*sqrt(2*pi)*Fnoise(2);
  Xem(j+1,:) = [x1 x2];
end

save("stochastic_new_periodic_ex2_1000_01.mat","t","Xem")

%% Visualization

load("stochastic_new_periodic_ex2_1000_01.mat","t","Xem")

for secn = [0.1 0.9]
  
  Y=[];
  c=0;
  for i=1:length(t)
    tmod = t(i)-floor(t(i));
    if abs(tmod-secn)<1e-8
      c = c+1;
      Y(c,:) = Xem(i,:); %#ok<SAGROW>
    end
  end
  
  % Visualization
  figure
  grid on
  box on
  hold on
  
  c2 = cos(2*pi*secn);
  s2 = sin(2*pi*secn);
  c4 = cos(4*pi*secn);
  s4 = sin(4*pi*secn);
  c8 = cos(8*pi*secn);
  s8 = sin(8*pi*secn);
  
  C = (1/32)*[4+c4-s4 -c4-s4;-c4-s4 4-3*c4-s4];
  
  lambda1 = (4-c4-s4-sqrt(3+2*c8+s8))/32;
  lambda2 = (4-c4-s4+sqrt(3+2*c8+s8))/32;
  
  [U,V] = eig(C);
  U1 = U(:,1); U2 = U(:,2);
  gamma = [s2; c2];
  
  tau = 0:2*pi/100:2*pi;
  
  for radius = 1:2
    x1 = zeros(1,length(tau));
    x2 = zeros(1,length(tau));
    for i=1:length(tau)
      x1(i) = radius*sigma*(U1(1)*sqrt(lambda1)*cos(tau(i)) ...
        +U2(1)*sqrt(lambda2)*sin(tau(i)))+gamma(1);
      x2(i) = radius*sigma*(U1(2)*sqrt(lambda1)*cos(tau(i)) ...
        +U2(2)*sqrt(lambda2)*sin(tau(i)))+gamma(2);
    end
    plot(x1, x2, 'k', 'LineWidth', 3)
  end
  plot(Y(:,1), Y(:,2), 'r.', 'Markersize',8)
  
  mu = mean(Y);
  cv = cov(Y);
  [W,D] = eig(cv);
  W1 = W(:,1); W2 = W(:,2);
  l1 = D(1,1); l2 = D(2,2);
  
  for radius = 1
    x1 = zeros(1,length(tau));
    x2 = zeros(1,length(tau));
    for i=1:length(tau)
      x1(i) = radius*(W1(1)*sqrt(l1)*cos(tau(i)) ...
        +W2(1)*sqrt(l2)*sin(tau(i)))+mu(1);
      x2(i) = radius*(W1(2)*sqrt(l1)*cos(tau(i)) ...
        +W2(2)*sqrt(l2)*sin(tau(i)))+mu(2);
    end
    plot(x1, x2, 'b--', 'LineWidth', 2)
  end
  plot(gamma(1), gamma(2), 'k.', 'MarkerFaceColor', 'k', 'MarkerSize', 12)
  
  set(gca, 'Linewidth', 2, 'Fontsize', 14)
  xlabel('$x_1$', 'Interpreter', 'Latex', 'Fontsize', 20)
  ylabel('$x_2$', 'Interpreter', 'Latex', 'Fontsize', 20)
  hold off
  if secn==0.9
    title('$t=0.9$', 'Interpreter', 'Latex', 'Fontsize', 18)
    axis([gamma(1)-.175 gamma(1)+.175 gamma(2)-0.175 gamma(2)+0.175])
  else
    title('$t=0.1$', 'Interpreter', 'Latex', 'Fontsize', 18)
    axis([gamma(1)-.15 gamma(1)+.15 gamma(2)-0.1 gamma(2)+0.1])
  end
end
