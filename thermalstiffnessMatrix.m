function [ke_t] = thermalstiffnessMatrix(a,b,t,kappa)

ke_t = kappa * t * ... 
[ a/(3*conj(b)) + b/(3*conj(a)),   a/(6*conj(b)) - b/(3*conj(a)), - a/(6*conj(b)) - b/(6*conj(a)),   b/(6*conj(a)) - a/(3*conj(b));
  a/(6*conj(b)) - b/(3*conj(a)),   a/(3*conj(b)) + b/(3*conj(a)),   b/(6*conj(a)) - a/(3*conj(b)), - a/(6*conj(b)) - b/(6*conj(a));
- a/(6*conj(b)) - b/(6*conj(a)),   b/(6*conj(a)) - a/(3*conj(b)),   a/(3*conj(b)) + b/(3*conj(a)),   a/(6*conj(b)) - b/(3*conj(a));
  b/(6*conj(a)) - a/(3*conj(b)), - a/(6*conj(b)) - b/(6*conj(a)),   a/(6*conj(b)) - b/(3*conj(a)),   a/(3*conj(b)) + b/(3*conj(a))];

end