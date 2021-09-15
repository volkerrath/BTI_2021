function [rc]=rcl(lambda)
% calculates volumetric heat capacity of rocks based on 
% a regression of diffusivity on thermal conductivity.
% (Ilmo Kukkonen, pers. comm.) 
A=0.5504;B=-0.1129;
rc=1.e6*lambda./(A*lambda+B);
