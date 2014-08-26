function answer=LocalToGlobal(nodes,localPt,Basis,varargin)
% answer=LocalToGlobal(nodes,localPt,Basis,varargin)

w=Basis(localPt(1),localPt(2),varargin{:});
answer=nodes(1,:)*w;