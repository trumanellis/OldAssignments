function elemEdgeMapping=ComputeElementEdgeMapping2(topo,EdgeList)

elemEdgeMapping=zeros(size(topo,1),size(topo,2));
temp=zeros(size(topo,1),1);
temp2=zeros(1,size(topo,1));
edgeTopo=[1,2;2,3;4,3;1,4];

for nE=1:size(topo,2)
%     temp=[];
    it=0;
    for nL=1:size(EdgeList,2)
        n1=EdgeList(1,nL);
        n2=EdgeList(2,nL);
        
        % Scan the edge list and add edge ID to element list as necessary
        if isempty(temp) && max(ismember(topo(:,nE),n1)) && max(ismember(topo(:,nE),n2))
            it=it+1;
            temp(it)=nL;
        elseif max(ismember(topo(:,nE),n1)) && max(ismember(topo(:,nE),n2))
            it=it+1;
            temp(it)=nL;
        end
    end
        
    % Now sort the list according to the local standard
%     temp2=[];
    it2=0;
    elem=topo(:,nE);
    oldData=EdgeList(:,temp(:));

    for iL=1:4
        edgeNodes=elem(edgeTopo(iL,:));
        [junk,iD]=intersect(oldData',edgeNodes','rows');
        it2=it2+1;
        temp2(it2)=iD;
    end
    sort=temp2(:);
    elemEdgeMapping(:,nE)=temp(sort);
end