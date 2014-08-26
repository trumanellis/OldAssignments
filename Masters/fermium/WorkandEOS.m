%% Work and EOS   
for i=1:NZ
    % Compute the new internal energy: F dot V internal energy update
    zonevel=0.5*(NEWvelocity(:,Quadmap(:,i))+OLDvelocity(:,Quadmap(:,i)));
    if HighOrderEnergy
        if StrongMass
            [cF]=WeightedCornerForce_J(allnodes(:,Quadmap(:,i)),QP,QuadWgts2d,...
                QuadGradVBasis,QuadPBasis,QuadInvJacobian(:,:,:,i),QuadDetJacobian(:,i),...
                ndofpZ,npdof,QuadPressure(:,i),QuadDensity(:,i),...
                OLDvelocity(:,Quadmap(:,i)),gamma,len(i),qlin,qquad,Qfrac);
        else
            [cF]=WeightedCornerForce(allnodes(:,Quadmap(:,i)),QP,QuadWgts2d,...
                QuadGradVBasis,QuadPBasis,QuadInvJacobian(:,:,:,i),QuadDetJacobian(:,i),...
                ndofpZ,npdof,QuadPressure(:,i),QuadDensity(:,i),...
                OLDvelocity(:,Quadmap(:,i)),gamma,len(i),qlin,qquad,Qfrac);
        end
        for j=1:npdof
            NEWenergy(j,i)=OLDenergy(j,i)-dt/cornerMass(j,i)*sum(dot(zonevel,cF(:,:,j)));
        end
    else
        NEWenergy(:,i)=OLDenergy(:,i)-dt/massZ(i)*sum(dot(zonevel,cornerForce(:,:,i)));
    end
    if ~StrongMass
		if npdof ~= 2
			% Compute new density
			localMP(:,:,i)=zeros(npdof,npdof);
			for n=1:QP
				localMP(:,:,i)=localMP(:,:,i)+QuadWgts2d(n)*refMP(:,:,n)*QuadDetJacobian(n,i);
			end
			NEWdensity(:,i)=localMP(:,:,i)\cornerMass(:,i);
			%EOS for an ideal monatomic gas
			pressure(:,i)=(gamma-1)*NEWdensity(:,i).*NEWenergy(:,i);
			% Check for negative pressure
			if min(pressure(:,i)) < -1e-6
				fprintf('WARNING -- Negative Pressure detected at cell %4.0f\n',i);
			elseif min(pressure(:,i)) < 0
				for s=1:npdof
					if pressure(s,i) < 0
						pressure(s,i) = 0;
					end
				end
			end
			% Evaluate density and pressure at quadrature points
			QuadDensity(:,i)=(NEWdensity(:,i)'*QuadPBasis)';
			QuadPressure(:,i)=(pressure(:,i)'*QuadPBasis)';
		else
			% Determine if even or odd zone
			if mod(i,2)
				ZoneIndex=1;
			else
				ZoneIndex=2;
			end
			% Compute new density
			localMP(:,:,i)=zeros(npdof,npdof);
			for n=1:QP
				localMP(:,:,i)=localMP(:,:,i)+QuadWgts2d(n)*refMP(:,:,n,ZoneIndex)*QuadDetJacobian(n,i);
			end
			NEWdensity(:,i)=localMP(:,:,i)\cornerMass(:,i);
			%EOS for an ideal monatomic gas
			pressure(:,i)=(gamma-1)*NEWdensity(:,i).*NEWenergy(:,i);
			% Check for negative pressure
			if min(pressure(:,i)) < -1e-6
				fprintf('WARNING -- Negative Pressure detected at cell %4.0f\n',i);
			elseif min(pressure(:,i)) < 0
				for s=1:npdof
					if pressure(s,i) < 0
						pressure(s,i) = 0;
					end
				end
			end
			% Evaluate density and pressure at quadrature points
			QuadDensity(:,i)=(NEWdensity(:,i)'*QuadPBasis(:,:,ZoneIndex))';
			QuadPressure(:,i)=(pressure(:,i)'*QuadPBasis(:,:,ZoneIndex))';
		end
    else
            % Constant Density*detJ
            %EOS for an ideal monatomic gas
%             pressure(:,i)=(gamma-1)*NEWdensity(:,i).*NEWenergy(:,i);
            pressure(:,i)=pressure(:,i).*NEWenergy(:,i)./OLDenergy(:,i);
            % Check for negative pressure
            if min(pressure(:,i)) < -1e-6
                fprintf('WARNING -- Negative Pressure detected at cell %4.0f\n',i);
            elseif min(pressure(:,i)) < 0
                for s=1:npdof
                    if pressure(s,i) < 0
                        pressure(s,i) = 0;
                    end
                end
            end
            % Evaluate pressure at quadrature points
            QuadPressure(:,i)=(pressure(:,i)'*QuadPBasis)';
    end
end
    
% Check energy conservation
if HighOrderEnergy
    internalenergy(cycle)=sum(cornerMass(:).*OLDenergy(:));
else
    internalenergy(cycle)=sum(massZ.*OLDenergy);
end
if ~FullMassMatrixSolve
    kineticenergy(cycle)=0.5*(OLDvelocity(1,:)*MMQuad*OLDvelocity(1,:)'+OLDvelocity(2,:)*MMQuad*OLDvelocity(2,:)');
    sumMassN(cycle)=sum(massN(:));
else
    kineticenergy(cycle)=0.5*(OLDvelocity(1,:)*MMX*OLDvelocity(1,:)'+OLDvelocity(2,:)*MMY*OLDvelocity(2,:)');
    sumMassN(cycle)=sum(MMQuad(:));
end
        
% Increment time
t=t+dt;
if t>=tstop
   stop=true;
end
timedata(cycle)=t;

% Reset old state variables to new state variables
OLDvelocity=NEWvelocity;
OLDdensity=NEWdensity;
OLDenergy=NEWenergy;