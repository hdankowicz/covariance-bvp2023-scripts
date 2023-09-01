function [out1,out2] = sym_bistable_rhs_0(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12)
%SYM_BISTABLE_RHS_0
%    [OUT1,OUT2] = SYM_BISTABLE_RHS_0(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    14-Aug-2023 10:51:31

out1 = in3;
if nargout > 1
    t2 = -in2;
    out2 = t2+in5.*cos((in1.*pi.*2.0)./in4)-in3.*in6+in2.^2.*t2;
end