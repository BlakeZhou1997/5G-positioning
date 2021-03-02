function [phi_out,f,fd]=Chapter_17_Function_6(v_hat,v_training,phi_hat_training);

f=[]; fd=[];
if nargin==1
  [phi_out,f,fd]=Chapter_17_Function_7(v_hat);
elseif nargin==3
  interp_method='spline';
  v_training2=[-v_training; v_training];
  phi_hat_training2=[-phi_hat_training; phi_hat_training];
  [v_training2 ind]=sort(v_training2);
 phi_out=interp1(v_training2,phi_hat_training2(ind),v_hat,interp_method,0);
  ind=find(phi_out==0);
  if ~isempty(ind)
    phi_out(ind)=abs(phi_hat_training2(find(v_training2== ...
        max(v_training2))))*sign(v_hat(ind));
  end;
end;

return;