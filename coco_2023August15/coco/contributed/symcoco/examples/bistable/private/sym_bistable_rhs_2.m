function [out1,out2] = sym_bistable_rhs_2(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12)
%SYM_BISTABLE_RHS_2
%    [OUT1,OUT2] = SYM_BISTABLE_RHS_2(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    14-Aug-2023 10:51:31

out1 = 0.0;
if nargout > 1
    t2 = 1.0./in4;
    t3 = t2.^2;
    t4 = in1.*t2.*pi.*2.0;
    t5 = in7.*t2.*pi.*2.0;
    t6 = in1.*in10.*t3.*pi.*2.0;
    t7 = sin(t4);
    t8 = -t6;
    t9 = t5+t8;
    out2 = in9.*in12.*-2.0-in2.*in8.^2.*6.0-in5.*t7.*(in1.*in10.^2.*t2.^3.*pi.*4.0-in7.*in10.*t3.*pi.*4.0)-in5.*t9.^2.*cos(t4)-in11.*t7.*t9.*2.0;
end
