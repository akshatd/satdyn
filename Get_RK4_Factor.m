function [K1,K2,K3,K4,K] = Get_RK4_Factor(f_xd,Xk,Uk,T,params)
% Xd = f_xd(X,U,params), column vectors. Compute K such that X(k+1) = X(k) + K(X(k),U(k))T

% Compute gradients
K1 = f_xd(Xk,Uk,params);
K2 = f_xd(Xk + K1*T/2,Uk,params);
K3 = f_xd(Xk + K2*T/2,Uk,params);
K4 = f_xd(Xk + K3*T,Uk,params);

K = (K1 + 2*K2 + 2*K3 + K4)/6;

end

