function refnodes=ComputeReferenceMeshNodes(NZx,NZy,xRange,yRange)

xmin=xRange(1);
xmax=xRange(2);
ymin=yRange(1);
ymax=yRange(2);

NN=(NZx+1)*(NZy+1);
refnodes=zeros(2,NN);
x=xmin:(xmax-xmin)/NZx:xmax;
y=ymin:(ymax-ymin)/NZy:ymax;

for i=1:(NZy+1)
    for j=1:(NZx+1)
        c=j+(i-1)*(NZx+1);
        refnodes(1,c)=x(j);
        refnodes(2,c)=y(i);
    end
end

