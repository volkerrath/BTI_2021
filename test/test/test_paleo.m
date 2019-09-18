% v.r. 8/01                    last change: Aug.16,2001 
year=31557600;
t=[0:1000:30000];nt=length(t);
t=t*year;

[T1,tn]=paleo_haenel(t);
figure;plot(tn/year,T1);
title('Haenels Paleoclimate');grid on;
xlabel('Time [a]');ylabel('Surface Temperature [\circ C]');


tt=10000*year;
[T2]=paleo_step(t,5,tt);
figure;plot(t/year,T2);
title('Step function');grid on;
xlabel('Time [a]');ylabel('Surface Temperature [\circ C]');

