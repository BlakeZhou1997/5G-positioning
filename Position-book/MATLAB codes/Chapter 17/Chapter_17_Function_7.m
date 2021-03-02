function [phi_akua,f_akus,fd_akua]=Chapter_17_Function_7(x,x_grid)

N=length(x);

kernel='gaussian';

% smoothing factor
alpha=0.5;

%x_grid=colvec(linspace(-max(abs(x)),max(abs(x)),M));
%delta_grid=mean(diff(x_grid));
if nargin==1
  x_grid=x;
end;
M=length(x_grid);

hR=0.79*iqr(x)*N^(-1/5);


%***** Adaptive Kernel Estimate *****% 
f_pilot=mean(feval([kernel,'_kernel'],(x*ones(1,N)-ones(N,1)*x.')/hR),2)/hR;
% Pilot estimate

lambda=(f_pilot/geomean(f_pilot)).^(-alpha);

f_ak=mean(feval([kernel,'_kernel'],(x_grid*ones(1,N)-ones(M,1)*x.')./(hR*ones(M,1)*lambda.'))./(hR*ones(M,1)*lambda.'),2);
%************************************%

  
%***** AKE with attempt to fix tails and ensure unimodality *****% 
% evalute at x and midpoints 
lambda_dash=fixtails(x,kernel,hR,lambda); 
% first attempt to fix the tails 
f_akt=mean(feval([kernel,'_kernel'],(x_grid*ones(1,N)-ones(M,1)*x.')./(hR*ones(M,1)*lambda_dash.'))./(hR*ones(M,1)*lambda_dash.'),2);
h=unimodal(x,kernel,hR,lambda_dash); % ensure unimodality
f_aku =mean(feval([kernel,'_kernel'],(x_grid*ones(1,N)-ones(M,1)*x.')./(h*ones(M,1)*lambda_dash.'))./(h*ones(M,1)*lambda_dash.'),2);
f_akum=mean(feval([kernel,'_kernel'],(x_grid*ones(1,N)+ones(M,1)*x.')./(h*ones(M,1)*lambda_dash.'))./(h*ones(M,1)*lambda_dash.'),2);
fd_aku = mean(feval([kernel,'_kernel_d'],(x_grid*ones(1,N)-ones(M,1)*x.')./(h*ones(M,1)*lambda_dash.'))./(h*ones(M,1)*lambda_dash.').^2,2);
fd_akum=-mean(feval([kernel,'_kernel_d'],(x_grid*ones(1,N)+ones(M,1)*x.')./(h*ones(M,1)*lambda_dash.'))./(h*ones(M,1)*lambda_dash.').^2,2);
f_akus =(f_aku+f_akum)/2;
fd_akua=(fd_aku-fd_akum)/2; 
%****************************************************************%

ind=find(f_akus==0);
fd_akua(ind)=0; f_akus(ind)=eps;
phi_akua=-fd_akua./f_akus;

return;
