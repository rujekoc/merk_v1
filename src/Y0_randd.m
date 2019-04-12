function Y = Y0_randd(x,ep)
% usage: Y = Y0_randd(x,ep)
% gives initial value for reaction diffusion test problem.
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% April 2019

lambda = .5*sqrt(2*ep*1e4);
Y = (1+ exp(lambda*(x-1))).^(-1);
Y = Y';

end
%end of script
