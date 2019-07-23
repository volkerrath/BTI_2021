classdef netcdfnamedlist < handle
    % Helper netcdfobj.m
    %
    %
    % Aslak Grinsted 2009
    properties(Hidden)
        listdata = [];
    end
    
    methods
        function obj = netcdfnamedlist(noitems)
            obj.listdata=cell(noitems,1);
        end
        
        function value=subsref(obj,s)
            for ii=1:length(obj.listdata)
                if strcmp(s(1).subs,obj.listdata{ii}.name)
                    value=obj.listdata{ii};
                    for jj=2:length(s)
                        value=subsref(value,s(jj));
                    end
                    return
                end
            end
            error(sprintf('''%s'' not found.',s.subs));           
        end
%         
%         function value=subsindex(obj,name)
%             for ii=1:length(listdata)
%                 if strcmp(name,listdata{ii}.name)
%                     value=ii;
%                     return
%                 end
%             end
%             error(sprintf('''%s'' not found.',name));
%         end
        
        
        function display(obj)
            for ii=1:length(obj.listdata)
                prettydisp(obj.listdata{ii});
            end
        end
        
    end
    

end
