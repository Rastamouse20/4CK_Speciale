function [A_e0] = CouplingMatrix(CTE,a,b,t,nu)
% returns the zero CouplingMatrix, A_e0, in the problem
%       f_thermal(t) = A_e * t_e
%       f_thermal(t) = E_e*A_e0 * t_e
% Where
%       f_thermal is the equvilent thermal load force
%       A_e is the coupling matrix
%       t_e is the temperature at the nodes
A_e0 = ...
[ (CTE*b*t)/(3*(nu - 1)),  (CTE*b*t)/(3*(nu - 1)),  (CTE*b*t)/(6*(nu - 1)),  (CTE*b*t)/(6*(nu - 1));
  (CTE*a*t)/(3*(nu - 1)),  (CTE*a*t)/(6*(nu - 1)),  (CTE*a*t)/(6*(nu - 1)),  (CTE*a*t)/(3*(nu - 1));
 -(CTE*b*t)/(3*(nu - 1)), -(CTE*b*t)/(3*(nu - 1)), -(CTE*b*t)/(6*(nu - 1)), -(CTE*b*t)/(6*(nu - 1));
  (CTE*a*t)/(6*(nu - 1)),  (CTE*a*t)/(3*(nu - 1)),  (CTE*a*t)/(3*(nu - 1)),  (CTE*a*t)/(6*(nu - 1));
 -(CTE*b*t)/(6*(nu - 1)), -(CTE*b*t)/(6*(nu - 1)), -(CTE*b*t)/(3*(nu - 1)), -(CTE*b*t)/(3*(nu - 1));
 -(CTE*a*t)/(6*(nu - 1)), -(CTE*a*t)/(3*(nu - 1)), -(CTE*a*t)/(3*(nu - 1)), -(CTE*a*t)/(6*(nu - 1));
  (CTE*b*t)/(6*(nu - 1)),  (CTE*b*t)/(6*(nu - 1)),  (CTE*b*t)/(3*(nu - 1)),  (CTE*b*t)/(3*(nu - 1));
 -(CTE*a*t)/(3*(nu - 1)), -(CTE*a*t)/(6*(nu - 1)), -(CTE*a*t)/(6*(nu - 1)), -(CTE*a*t)/(3*(nu - 1))];
end