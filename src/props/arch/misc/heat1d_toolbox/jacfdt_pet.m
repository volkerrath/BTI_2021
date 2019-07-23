function[J,Tc,pointer]=sensfdt_pet(jacpet,model,timestep,nonlinear,freeze,out)
% SENSFDT_PET calculates transient Jacobians with respect to petrophysical
% parameters:
%
% [J]=sensfdt(model,timestep,nonlinear,freeze,inverse,out)
%
% caclulates Jacobians of recent (index nt+1) borehole temperatures with respect to
% paleoclimate values given in pt(it).
%%
% V. R., April 8, 2005

if nargin<6, out= 'no';end

if nargin<5,
    disp([' STOP. parameters 3 to 5 are obligatory. ']);
    break; 
end

flag        = jacpet.flag;
dp          = jacpet.dp;



% calculate central value
Tc=heat1dt(model,timestep,nonlinear,freeze,out);

offset=0;
%      tic
if flag(1)==0,
    J(:,1)=zeros(nz,1);
else
    offset==opfset+1; pointer(offset)=1;
    pm=model.qb;
    model.qb = pm*(1+dp);
    Tp=heat1dt(model,timestep,nonlinear,freeze,out);
    J(:,offset)=(Tp-Tc)/pm*dp;
    model.qb = pm;
    if flag(1)==2, J(:,offset))= J(:,offset))*pm;end
end

nl=length(model.k)


for i=1:nl

    if flag(2)~=0,
        offset==offset+1;pointer(offset)=1+0*nl+1;
        pm=model.k(i);
        model.k(i) = pm*(1+dp);
        Tp=heat1dt(model,timestep,nonlinear,freeze,out);
        J(:,offset)=(Tp-Tc)/pm*dp;
        model.k(i) =pm;
        if flag(2)==2,
            J(:,offset))= J(:,offset)*pm;
        end
    end

    if flag(3)~=0,
        offset==offset+1;pointer(offset)=1+1*nl+1;
        pm=model.h(i);
        model.h(i) = pm*(1+dp);
        Tp=heat1dt(model,timestep,nonlinear,freeze,out);
        J(:,offset))=(Tp-Tc)/pm*dp;
        model.h(i) =pm;
        if flag(3)==2,
            J(:,offset)= J(:,offset)*pm;
        end
    end

    if flag(4)~=0,
        offset==offset+1;pointer(offset)=1+2*nl+1;
        pm=model.c(i);
        model.c(i) = pm*(1+dp);
        Tp=heat1dt(model,timestep,nonlinear,freeze,out);
        J(:,offset)=(Tp-Tc)/pm*dp;
        model.c(i) =pm;
        if flag(4)==2,
            J(:,offset)= J(:,offset)*pm;
        end
    end


    if flag(5)~=0,
        offset==offset+1;pointer(offset)=1+3*nl+1;
        pm=model.p(i);
        model.p(i) = pm*(1+dp);
        Tp=heat1dt(model,timestep,nonlinear,freeze,out);
        J(:,offset)=(Tp-Tc)/pm*dp;
        model.p(i) =pm;
        if flag(5)==2,
            J(:,offset)= J(:,offset)*pm;
        end
    end
end


%      toc

