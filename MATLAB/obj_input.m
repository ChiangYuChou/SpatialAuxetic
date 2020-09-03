function [v l f] = obj_input(filename)
%OBJ_INPUT Loading obj to v:ver l:indexs of line f:indexs of face
%   example:
%	object=fopen('obj.obj','r');
%	[v l f]=obj_input(object);
%	obj_draw(v,l)
object=fopen(filename,'r');
iv=0;
il=0;
ifa=0;
v=[0 0 0];
l=[0 0];
f={[0 0 0]};
while ~feof(object)
    line=fgetl(object);
    if isempty(line)
    elseif line(1)=='v'
        if line(2)==' '
            iv=iv+1;
            v(iv,:)=sscanf(line(2:end),'%f',3);
        end
    elseif line(1)=='l'
        if line(2)==' '
            il=il+1;
            l(il,:)=sscanf(line(2:end),'%d',2);
        end
    elseif line(1)=='f'
        if line(2)==' '
            ifa=ifa+1;
            fline=findstr(line, ' ');
            for i=1:length(fline)
                j=sscanf(line(fline(i):end),'%d',1);
                if j
                f{ifa}(i)=j;
                end
            end
        end
    end
end

end

