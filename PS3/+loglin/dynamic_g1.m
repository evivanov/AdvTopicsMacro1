function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = loglin.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(3, 8);
g1(1,3)=(-1)/(y(3)*y(3));
g1(1,6)=(-(params(2)*(-T(1))/(y(6)*y(6))));
g1(1,4)=(-(params(2)*exp(y(7))*params(1)*getPowerDeriv(y(4),params(1)-1,1)/y(6)));
g1(1,7)=(-(params(2)*T(1)/y(6)));
g1(2,3)=1;
g1(2,1)=(-(exp(y(5))*getPowerDeriv(y(1),params(1),1)));
g1(2,4)=1;
g1(2,5)=(-T(2));
g1(3,2)=(-(params(3)*1/y(2)));
g1(3,5)=1/y(5);
g1(3,8)=(-params(4));

end
