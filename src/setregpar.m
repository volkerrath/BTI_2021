function [regpar]= setregpar(regvec0,reg0vec,reg1vec,reg2vec)
% generates matrix of regularization vecamters regvec
            ireg=0;
            for ireg0=1:length(reg0vec)
                for ireg1=1:length(reg1vec)
                    for ireg2=1:length(reg2vec)
                        ireg=ireg+1;
                        regpar(ireg,1)=reg0vec(ireg0)*regvec0(1);
                        regpar(ireg,2)=reg1vec(ireg1)*regvec0(2);
                        regpar(ireg,3)=reg2vec(ireg2)*regvec0(3);
                    end
                end
            end
end