function xArr = ODE_RK4(f_xd,x0,N,uArr,T,params)
% Solve IVP for nonlinear dynamics, discretize using RK4
%  ---------- Inputs -------------
% f_xd = @(X,U,params)[eye(2)*x + u]; % Function handle for nonlinear dynamics xdot = f_xd(x,u) as column vector
% x0 = [1;2]; % Initial state column vector
% N = 10; % Number of steps in the prediction horizon
% uArr = ones(2,N); % Control trajectory for u0 ... u(N-1), where uk is a column vector
% T = 0.1; % Sampling period

%  -------------- Outputs -----------------
% xArr = column vector of x0 (given), x1, .., xN

%% Obtain dimention of state and inputs
ns = length(x0);
nu = size(uArr,1);

%% Allocate memory
xArr = zeros(ns,N+1);
xArr(:,1) = x0;

%% RK4 
for n=0:(N-1)    
    % Get current state and input
   xk = xArr(:,n+1); uk = uArr(:,n + 1);
   
   [~,~,~,~,K] = Get_RK4_Factor(f_xd,xk,uk,T,params);
   % Propagate
   xArr(:,n + 1 + 1) = xArr(:,n  + 1) + K*T;
end