function [x] = cauchy(a,b,d)
%CAUCHY �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
  %a-x0=0.b-gamma=0.5...��任�����������ֲ���C��x0,gamma��
   u=rand(1,d);
   x=a-b./tan(3.1415926.*u);
end

