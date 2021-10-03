function [y] = f4(x1,x2)
%laser cutting
y=zeros(size(x1));

yval=0.01*( ( (x1-7).^2+(x2-7).^2-62 ).^2+( (x1-7).^2-0.5*(x1-7)-0.5*(x2-7)-1.5 ).^2 )+10;
y(x1.^2+x2.^2>80)=yval(x1.^2+x2.^2>80);
end

