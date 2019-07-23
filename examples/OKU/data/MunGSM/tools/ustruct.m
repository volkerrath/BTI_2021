function [outstruct]=ustruct(ostruct,nstruct)
% uses the structure nstruct to add or overwrite fields of original
% structure ostruct.
% vr 7/15 (based on overwrite_merge by Damial Amir, 2012)
if (isstruct(nstruct)&& isstruct(ostruct))
    fieldso=fields(nstruct);length_o=numel(fieldso);
    fieldsn=fields(ostruct);length_o=numel(fieldso);
    
    for fieldo=1:length_o
        fieldnameo= fieldso{fieldo};
        found=0;
        for fieldn=1:length_n
            fieldnamen=fieldsn{fieldn};
            if (strcmp(fieldnameo,fieldnamen))
                found=1;
                ostruct.(fieldnamen)=nstruct.(fieldnameo);
            end 
        end 
        
        if (found==0 )
            ostruct.(fieldnameo)=nstruct.(fieldnameo);
            
        end
        
    end 
end 

outstruct=ostruct;


end