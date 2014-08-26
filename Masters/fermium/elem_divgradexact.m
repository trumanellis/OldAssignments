function cF=elem_divgradexact(nodes,zonevel,QP,wgts2d,pts2d,VBasis,GradVBasis,invJ,detJ,ndofpZ,m,n)

% % Compute corner forces
cF=zeros(2,ndofpZ);
% 
% for n=1:QP
%     % Grad p Corner Force
%     globx=LocalToGlobal(nodes(1,:),pts2d(n,:),VBasis);
%     globy=LocalToGlobal(nodes(2,:),pts2d(n,:),VBasis);
% 	gradVx=[-m*pi*sin(m*pi*globx)*cos(n*pi*globy), -n*pi*cos(m*pi*globx)*sin(n*pi*globy)];
%     gradVy=[-n*pi*sin(n*pi*globx)*cos(m*pi*globy), -m*pi*cos(n*pi*globx)*sin(m*pi*globy)];
% %     gradVx=[pi/2*cos(pi/2*globx), 0];
% %     gradVy=[0, pi/2*cos(pi/2*globy)];
%     stiffmat=stiffmat-wgts2d(n)*(invJ(:,:,n)*[gradVx;gradVy])*(invJ(:,:,n)*[gradVx;gradVy])'*detJ(n);
% end
% cF=-stiffmat*zonevel;