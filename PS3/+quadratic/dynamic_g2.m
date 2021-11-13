function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
% function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
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
%   g2
%

if T_flag
    T = quadratic.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
v2 = zeros(16,3);
v2(1,1)=1;
v2(2,1)=1;
v2(3,1)=1;
v2(4,1)=1;
v2(5,1)=1;
v2(6,1)=1;
v2(7,1)=1;
v2(8,1)=1;
v2(9,1)=1;
v2(10,1)=1;
v2(11,1)=2;
v2(12,1)=2;
v2(13,1)=2;
v2(14,1)=2;
v2(15,1)=3;
v2(16,1)=3;
v2(1,2)=19;
v2(2,2)=46;
v2(3,2)=44;
v2(4,2)=30;
v2(5,2)=47;
v2(6,2)=54;
v2(7,2)=28;
v2(8,2)=31;
v2(9,2)=52;
v2(10,2)=55;
v2(11,2)=1;
v2(12,2)=5;
v2(13,2)=33;
v2(14,2)=37;
v2(15,2)=10;
v2(16,2)=37;
v2(1,3)=(y(3)+y(3))/(y(3)*y(3)*y(3)*y(3));
v2(2,3)=(-(params(2)*(-((-T(1))*(y(6)+y(6))))/(y(6)*y(6)*y(6)*y(6))));
v2(3,3)=(-(params(2)*(-T(4))/(y(6)*y(6))));
v2(4,3)=v2(3,3);
v2(5,3)=(-(params(2)*(-T(1))/(y(6)*y(6))));
v2(6,3)=v2(5,3);
v2(7,3)=(-(params(2)*exp(y(7))*params(1)*getPowerDeriv(y(4),params(1)-1,2)/y(6)));
v2(8,3)=(-(params(2)*T(4)/y(6)));
v2(9,3)=v2(8,3);
v2(10,3)=(-(params(2)*T(1)/y(6)));
v2(11,3)=(-(exp(y(5))*getPowerDeriv(y(1),params(1),2)));
v2(12,3)=T(3);
v2(13,3)=v2(12,3);
v2(14,3)=(-T(2));
v2(15,3)=(-(params(3)*(-1)/(y(2)*y(2))));
v2(16,3)=(-1)/(y(5)*y(5));
g2 = sparse(v2(:,1),v2(:,2),v2(:,3),3,64);
end
