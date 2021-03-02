function qinv = qinverse(q)
% qinv = qinverse(q)
% Output the inverse of quarternion vector
% q = [q_vector; q_scalar]

qinv = [-q(1:3);q(4)];