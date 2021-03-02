% Compute pairwise potential - y(x1,x2):
% d12 - distance between x1 and x2, sigma - std. of distance measurements
function y=pairwise_pot(x1,x2,d12,R,sigma)
epsilon=0.000001; % some very small number
logy=0.5*(abs(x2-x1)-d12).^2/sigma^2; % log of the potential
min_logy=min(logy);
y=exp(-(logy-min_logy))+epsilon; % compute pairwise potential (and avoid y=0)
y=y/sum(y); % normalize
end
