function y=cpm(x)
% Y=CPM(X), cross product matrix
% input either a 3x1 matrix 1x3 matrix, will return a 3x3 cross product
% matrix of that particular matrix
%


y=[0 -x(3) x(2);
   x(3) 0 -x(1);
   -x(2) x(1) 0];



