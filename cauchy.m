function [x] = cauchy(a,b,d)
%CAUCHY 此处显示有关此函数的摘要
%   此处显示详细说明
  %a-x0=0.b-gamma=0.5...逆变换法产生柯西分布，C（x0,gamma）
   u=rand(1,d);
   x=a-b./tan(3.1415926.*u);
end

