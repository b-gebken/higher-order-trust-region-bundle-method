% A script for visualizing the first two derivatives of f for different p
% (for numerically verifying the stated properties)

x_sym = sym('x'); assume(x_sym,'real');

p = 4;

f_sym = x_sym^(p+1) * sin(1/x_sym) + 1/p * abs(x_sym)^p;

d1_f_sym = diff(f_sym,x_sym,1);
d2_f_sym = diff(f_sym,x_sym,2);

f = matlabFunction(f_sym,'Vars',{x_sym(:)});
d1_f = matlabFunction(d1_f_sym,'Vars',{x_sym(:)});
d2_f = matlabFunction(d2_f_sym,'Vars',{x_sym(:)});

x_arr = linspace(-0.1,0.1,10000);
tmp = [-1./(2*pi*(2:10000)),1./(2*pi*(2:10000))];
x_arr = sort([x_arr,tmp]);

figure
subplot(1,3,1)

plot(x_arr,f(x_arr),'k-');
hold on
yline(0,'r-')

axis square

subplot(1,3,2)

plot(x_arr,d1_f(x_arr),'k-');
hold on
yline(0,'r-')

axis square

subplot(1,3,3)

plot(x_arr,d2_f(x_arr),'k-');
hold on
yline(0,'r-')

axis square