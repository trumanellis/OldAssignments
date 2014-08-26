function edgeList=ComputeEdgeList(topo)

edgeList=[];
edgeTopo=[1,2;2,3;4,3;1,4];

for nE=1:size(topo,2)
    elem=topo(:,nE);
    for i=1:4
        edgeNodes=elem(edgeTopo(i,1:2));
        if  isempty(edgeList)
            edgeList=[edgeList,edgeNodes];
        elseif max(~ismember(edgeList,edgeNodes))
            edgeList=[edgeList,edgeNodes];
        end
    end
end

