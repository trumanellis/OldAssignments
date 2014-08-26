function centroid=ComputeCentroid(nodes)
% centroid=ComputeCentroid(nodes)
centroid=[sum(nodes(1,:));sum(nodes(2,:))]/length(nodes);