function [out1,out2] = sym_bistable_rhs_3(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12)
%SYM_BISTABLE_RHS_3
%    [OUT1,OUT2] = SYM_BISTABLE_RHS_3(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    14-Aug-2023 10:51:31

out1 = 0.0;
if nargout > 1
    t2 = in10.^2;
    t3 = 1.0./in4;
    t4 = t3.^2;
    t5 = t3.^3;
    t6 = in1.*t3.*pi.*2.0;
    t7 = in7.*t3.*pi.*2.0;
    t8 = cos(t6);
    t9 = in1.*in10.*t4.*pi.*2.0;
    t10 = in7.*in10.*t4.*pi.*4.0;
    t11 = sin(t6);
    t13 = in1.*t2.*t5.*pi.*4.0;
    t12 = -t9;
    t14 = -t13;
    t15 = t7+t12;
    t16 = t10+t14;
    out2 = in8.^3.*-6.0+in5.*t11.*t15.^3-in11.*t8.*t15.^2.*3.0-in5.*t11.*(in7.*t2.*t5.*pi.*1.2e+1-in1.*in10.^3.*t4.^2.*pi.*1.2e+1)+in11.*t11.*t16.*3.0+in5.*t8.*t15.*t16.*3.0;
end
