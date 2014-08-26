function answer=ComputeInvLengthSquared(nodes)

x=nodes(1,:);
y=nodes(2,:);

area=ComputeVolume(nodes);
dx=x(1)+x(2)-x(3)-x(4);
dy=y(1)+y(2)-y(3)-y(4);
dsqi=0.25*(dx*dx+dy*dy);
dx=x(2)+x(3)-x(4)-x(1);
dy=y(2)+y(3)-y(4)-y(1);
dsqj=0.25*(dx*dx+dy*dy);

answer=(dsqi+dsqj)/(area*area);