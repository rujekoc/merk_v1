function [A] = A_randd(a,b,n,ep)
% Usage: [A] = A_randd(a,b,n,ep)
% INPUTS : n - number of intervals in discretization
%          a - left end-point
%          b - right end-point
%         ep - epsilon value
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% April 2019

% set spatial mesh size
h = (b-a)/n;

% create matrix storage
A = spalloc(n+1,n+1,3*n+3);

% interior rows
for i=2:n
	A(i,i-1) = 1/h/h;
	A(i,i) = -2/h/h;
	A(i,i+1) = 1/h/h;
end

% boundary rows
A(1,1:2) = [-2 2]/h^2;
A(n+1,n:n+1) = [2 -2]/h^2;

A = 1e-4/ep*A;
end
% end of function
