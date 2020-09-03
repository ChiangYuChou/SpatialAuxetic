function [ ] = obj_draw(v,l,fig,config)
%OBJ_DRAW Summary of this function goes here
%   Detailed explanation goes here
figure(fig)
plot3(v(l(1,:),1),v(l(1,:),2),v(l(1,:),3),config); hold on
for i=2:length(l(:,1))
	plot3(v(l(i,:),1),v(l(i,:),2),v(l(i,:),3),config)
end
hold off;
axis equal

end

