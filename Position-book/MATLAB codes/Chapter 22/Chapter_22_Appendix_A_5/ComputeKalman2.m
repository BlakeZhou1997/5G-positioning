function outputs = ComputeKalman2(xn_1,Pn_1,Q_n1,H,R,Z)

P_stateOnly = Pn_1; % JTM 3/3/11
% Kalman gain
K = P_stateOnly * (H.')*inv(H * P_stateOnly * (H.') + R );
% Error covariance update
P = (eye(length(xn_1)) - K * H)* P_stateOnly;
% State update
X = K*Z; % see notes 1/2/11
%% Residuals
outputs.res = Z;

outputs.X = X;
outputs.P = P;

end