function [cpf,enthalp]=cpfT(T)
% nach Phillips 1981
%      entab in this case contains enthalpies in J/Kg for a 0.5 molal NACL
%      solution from zero to 300 degrees C at increments of 25 C.
      delT=25.;
      entab=[0.,24.044,  48.184,  72.375,  96.623, ...
             121.111, 145.800, 170.821, 196.289, 222.298,...
                 249.024, 276.674, 305.341]*4187.; 
%     Set enthalpy and heat capacity of fluid.
      [n1 n2]=size(T);if n1==1, T=T'; end 
      isub=floor(T/delT)+1;

      cpf=(entab(isub+1)-entab(isub))/delT;
      if nargout>1,
        ethalp=entab(isub)+(T-isub*delT)*cpf;
      end
