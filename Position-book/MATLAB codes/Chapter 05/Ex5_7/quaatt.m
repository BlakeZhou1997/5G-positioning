function C=quaatt(q)



Phi = [q(4)*eye(3)+cpm(q(1:3)); -q(1:3)'];
Psi = [q(4)*eye(3)-cpm(q(1:3)); -q(1:3)'];

C=Phi'*Psi;