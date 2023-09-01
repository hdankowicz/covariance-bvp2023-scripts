function varargout=dhopf(action,varargin)
%% Automatically generated with matlabFunction
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'nargs'
   varargout{1}=2;
   return
  case 'nout'
   varargout{1}=2;
   return
  case 'argrange'
   varargout{1}=struct('x',1:2,'p',3:6);
   return
  case 'argsize'
   varargout{1}=struct('x',2,'p',4);
   return
  case 'vector'
   varargout{1}=struct('x',1,'p',1);
   return
  case 'extension'
   varargout{1}='rhs';
   return
  case 'maxorder'
   varargout{1}=3;
   return
end
nout=2;
order=varargin{1};
f=str2func(sprintf('dhopf_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});
end



function [out1,out2] = dhopf_rhs_0(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12)
%DHOPF_RHS_0
%    [OUT1,OUT2] = DHOPF_RHS_0(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    11-Aug-2023 09:58:22

t2 = in1.^2;
t3 = in2.^2;
t4 = t2+t3;
out1 = in1.*in3-in2.*in4+t4.*(in1.*in5-in2.*in6);
if nargout > 1
    out2 = in1.*in4+in2.*in3+t4.*(in1.*in6+in2.*in5);
end
end


function [out1,out2] = dhopf_rhs_1(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12)
%DHOPF_RHS_1
%    [OUT1,OUT2] = DHOPF_RHS_1(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    11-Aug-2023 09:58:22

t2 = in1.^2;
t3 = in2.^2;
t4 = in1.*in7.*2.0;
t5 = in2.*in8.*2.0;
t6 = t2+t3;
t7 = t4+t5;
out1 = in1.*in9+in3.*in7-in2.*in10-in4.*in8+t7.*(in1.*in5-in2.*in6)+t6.*(in1.*in11+in5.*in7-in2.*in12-in6.*in8);
if nargout > 1
    out2 = in1.*in10+in2.*in9+in3.*in8+in4.*in7+t7.*(in1.*in6+in2.*in5)+t6.*(in1.*in12+in2.*in11+in5.*in8+in6.*in7);
end
end


function [out1,out2] = dhopf_rhs_2(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12)
%DHOPF_RHS_2
%    [OUT1,OUT2] = DHOPF_RHS_2(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    11-Aug-2023 09:58:22

t2 = in1.^2;
t3 = in2.^2;
t4 = in7.^2;
t5 = in8.^2;
t6 = in1.*in7.*2.0;
t7 = in2.*in8.*2.0;
t8 = t4.*2.0;
t9 = t5.*2.0;
t10 = t2+t3;
t11 = t6+t7;
t12 = t8+t9;
out1 = in7.*in9.*2.0-in8.*in10.*2.0+t12.*(in1.*in5-in2.*in6)+t10.*(in7.*in11.*2.0-in8.*in12.*2.0)+t11.*(in1.*in11+in5.*in7-in2.*in12-in6.*in8).*2.0;
if nargout > 1
    out2 = in7.*in10.*2.0+in8.*in9.*2.0+t12.*(in1.*in6+in2.*in5)+t10.*(in7.*in12.*2.0+in8.*in11.*2.0)+t11.*(in1.*in12+in2.*in11+in5.*in8+in6.*in7).*2.0;
end
end


function [out1,out2] = dhopf_rhs_3(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12)
%DHOPF_RHS_3
%    [OUT1,OUT2] = DHOPF_RHS_3(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    11-Aug-2023 09:58:22

t2 = in7.^2;
t3 = in8.^2;
t4 = in1.*in7.*2.0;
t5 = in2.*in8.*2.0;
t6 = t2.*2.0;
t7 = t3.*2.0;
t8 = t4+t5;
t9 = t6+t7;
out1 = t8.*(in7.*in11.*2.0-in8.*in12.*2.0).*3.0+t9.*(in1.*in11+in5.*in7-in2.*in12-in6.*in8).*3.0;
if nargout > 1
    out2 = t8.*(in7.*in12.*2.0+in8.*in11.*2.0).*3.0+t9.*(in1.*in12+in2.*in11+in5.*in8+in6.*in7).*3.0;
end
end
