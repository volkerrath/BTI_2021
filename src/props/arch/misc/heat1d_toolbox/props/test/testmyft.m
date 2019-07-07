clear all;close all;clc
T=[1:100];
[M1 T]=myfT1(T);
plot(T,M1,'-r');hold on

[M2 T]=myfT2(T);
plot(T,M2,':b');hold on

dM=M1-M2;