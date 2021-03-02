function qout = qmultiply(q1,q2)
% qout = q1 x q2
% apply for all quarternion multiplication
% q1 = [q1_vector;q1_scalar]
% q2 = [q2_vector;q2_scalar]
%
% Require matlab file
%     cpm.m

Psi = [q1(4)*eye(3)-cpm(q1(1:3)); -q1(1:3)'];
qout = [Psi q1]*q2;
