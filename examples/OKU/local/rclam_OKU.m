function [rc]=rcl(lambda,plot)
% calculates volumetric heat capacity of rocks based on 
% a regression of diffusivity on themal conuctivity.
% (Ilmo kukkonen, pers. comm.) 
A=0.503;B=0.0839;
rc=1.e6*lambda./(A*lambda+B)
