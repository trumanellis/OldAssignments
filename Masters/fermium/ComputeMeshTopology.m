function topo=ComputeMeshTopology(NZx,NZy)

topo=zeros(4,NZx*NZy);

for j=0:NZy-1
    for i=0:NZx-1
        topo(:,(i+1)+j*NZx)=[
            j*(NZx + 1) + 1 + i,
            j*(NZx + 1) + 2 + i,
            (j + 1)*(NZx + 1) + 2 + i,
            (j + 1)*(NZx + 1) + 1 + i];
disp(' ');
    end
end
