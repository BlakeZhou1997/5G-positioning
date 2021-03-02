function v = CRLB(X, sigma2)
% X the first row is the unknown-position sensor,the 2:end row are the
% anchors
k = size(X, 1);
M = k - 1;
d = sqrt(sum((X(2:end, :) - ones(M, 1)*X(1, :)).^2,2));
N = size(X,2);
D = (ones(M, 1)*X(1, :) - X(2:end, :))./(d*ones(1, N));
FIM = D'*diag(1./sigma2)*D;
H = inv(FIM);
v = trace(H);
