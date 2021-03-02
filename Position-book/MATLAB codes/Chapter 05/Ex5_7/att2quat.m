function qf = att2quat(Q)
%This will give the quarternion of Q matrix







q4 = sqrt(1+trace(Q))/2;
% 
if abs(q4) > 0
    q1 = 1/4/q4*(Q(2,3)-Q(3,2));
    q2 = 1/4/q4*(Q(3,1)-Q(1,3));
    q3 = 1/4/q4*(Q(1,2)-Q(2,1));
    qf = [q1;q2;q3;q4];
else
    Q = Q';
    vx = [1 0 0]';
    vy = [0 1 0]';
    ax1 = cross(vx,Q(:,1));
    ax1 = ax1/norm(ax1);
    ang = acos(Q(1,1));
    q1 = [ax1*sin(ang/2);cos(ang/2)];
    vyr = quaatt(q1)'*vy;
    vyr  = vyr/norm(vyr);
    % qq = [psifun(q1) q1]*[vy;0];
    % vyr = [psifun(qq) qq]*[-q1(1:3);q1(4)]
    % vyr = [psifun(qq) qq]*[q1(1:3);-q1(4)]
    % vyr = vyr(1:3);
    ang2 = acos(vyr'*Q(:,2));
    tvec = cross(vyr,Q(:,2));
    tvec = tvec/norm(tvec);
    tval = tvec'*Q(:,1)/norm(Q(:,1));
    q2 = [sign(tval)*[1 0 0]'*sin(ang2/2);cos(ang2/2)];
    q2 = q2/norm(q2);
    
    qf = [psifun(q2) q2]*q1;
end
