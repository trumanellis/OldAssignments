function area=ComputeVolume(nodes)

if (nargin <2)
    planar=1;
end

x=nodes(1,:);
y=nodes(2,:);

A124=0.5*((x(2)-x(1))*(y(4)-y(1))-(x(4)-x(1))*(y(2)-y(1)));
A234=0.5*((x(4)-x(3))*(y(2)-y(3))-(x(2)-x(3))*(y(4)-y(3)));

if planar==1
    area=A124+A234;
else
    area=2*pi*(A124*1/3*(y(1)+y(2)+y(4))+A234*1/3*(y(2)+y(3)+y(4)));
end
