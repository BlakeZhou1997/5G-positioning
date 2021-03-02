function [theta_hat_final,theta_hat,fv_final]=Chapter_17_Function_5(x,S,...
    theta0,v_training,max_iters,epsilon)
% Inputs
%  
% x : observations
% S : code matrix
% theta0 : initial estimate of theta
% v_training : training sequence (optional, default=[])
% max_iters : maximum number of iterations (optional, default=50) 
% epsilon : convergence criterion (optional, default=1e-4) %  
% Outputs
%
% theta_hat_final : final estimate of theta
% theta_hata : matrix, the columns being successive estimates of theta

global counter_rd; 
  
max_iters_DEFAULT=50;
epsilon_DEFAULT=1e-2;

if ~exist('v_training')
  training=0;
elseif isempty(v_training)
  training=0;
else
  training=1;
end;

if ~exist('max_iters')
  max_iters=max_iters_DEFAULT;
elseif isempty(max_iters)
  max_iters=max_iters_DEFAULT;
end;

if ~exist('epsilon')
  epsilon=epsilon_DEFAULT;
elseif isempty(epsilon)
  epsilon=epsilon_DEFAULT;
end;

% MP pseudo-inverse of the code matrix
theta_hat(:,1)=theta0(:); % initial estimate of theta S_pind=inv(S'*S)*S';

if training
  phi_hat_training=Chapter_17_Function_6(v_training);
end;

i=1; v_hat_org=x-S*theta_hat(:,i);
while i<max_iters

  % Find the residuals
  v_hat=x-S*theta_hat(:,i);
  
  % Estimate the score function
  if training
    phi_hat=Chapter_17_Function_6(v_hat,v_training,phi_hat_training);
  else
    phi_hat=Chapter_17_Function_6(v_hat);
  end;
  
  [ind ind]=sort(v_hat);
  
  mu(i)=1.25*max(abs(diff(phi_hat(ind))./diff(v_hat(ind))));
  %1/mu(i)*inv(S'*S)*S'*phi_hat
  theta_hat(:,i+1)=theta_hat(:,i)+1/mu(i)*inv(S'*S)*S'*phi_hat;
  
  
  if abs((theta_hat(:,i+1)-theta_hat(:,i))./(theta_hat(:,i)))/ ...
          length(theta_hat(:,1))<epsilon
    break;
  end;
  
  i=i+1;

end;

%figure; plot(v_hat_org,phi_hat); figure; plot(1:50,phi_hat);  %AA = [v_hat_org phi_hat.']
%%phi_hat_sort = sort(phi_hat); figure; plot(-25:24,phi_hat_sort);

theta_hat_final=theta_hat(:,end);

v_hat_final=x-S*theta_hat_final;
[phi_hat_final,fv_final,dfv_final]=Chapter_17_Function_6(v_hat_final);

counter_rd=[counter_rd i];

return;

