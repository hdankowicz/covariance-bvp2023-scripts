function [out1,out2] = sym_bistable_rhs_1(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12)
%SYM_BISTABLE_RHS_1
%    [OUT1,OUT2] = SYM_BISTABLE_RHS_1(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    14-Aug-2023 10:51:31

out1 = in9;
if nargout > 1
    t2 = 1.0./in4;
    t3 = in1.*t2.*pi.*2.0;
    out2 = -in8-in3.*in12-in6.*in9-in2.^2.*in8.*3.0+in11.*cos(t3)-in5.*sin(t3).*(in7.*t2.*pi.*2.0-in1.*in10.*t2.^2.*pi.*2.0);
end
