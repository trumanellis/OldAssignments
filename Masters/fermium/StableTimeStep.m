%% Stable Time Step
% Compute the length scale
len=ones(1,NZ);
for N=1:NZ
    for i=1:ndofpZ
        xdist=allnodes(1,Quadmap(i,N))-allnodes(1,Quadmap(i+1:ndofpZ,N));
        ydist=allnodes(2,Quadmap(i,N))-allnodes(2,Quadmap(i+1:ndofpZ,N));
        len(N)=min([len(N),sqrt(xdist.^2+ydist.^2)]);
    end
%     len(N)=1/sqrt(ComputeInvLengthSquared(oldnodes(:,topo(:,N))));
end

% Compute stable time increment
if cycle > 1
    % Ramp up comparison time so that it eventually stops being chosen
    dt=dtInit*cycle^1.2;
    dt=min([dt,.25*min(len)/max(sqrt(OLDvelocity(1,:).^2+OLDvelocity(2,:).^2)),...
        tplot(pcounter)-t,tsave(scounter)-t,dtMax]);
%     for i=1:NZ
%         zonevel=OLDvelocityQuad(:,Quadmap(:,i));
%         % Maximum velocity in zone
%         umax=sqrt(max(abs(zonevel(1,:)))^2+max(abs(zonevel(2,:)))^2);
%         dt=min([dt,0.25*len(i)/(umax+1e-6),dtMax]);
%     end
%     dt=min([dt,tplot(pcounter)-t,tsave(scounter)-t]);
else
    dt=dtInit;
    if dt < 1e-9
        stop = true;
    end
end

% Compute sound speed
% soundSpeedZ=sqrt(gamma*pressureZ./OLDdensityZ);