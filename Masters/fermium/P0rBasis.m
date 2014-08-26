function [w]=P0rBasis(x,y,index)

if mod(index,2)
    if y == x
        w=[.5;.5];
    elseif y < x
        w=[1;0];
    else
        w=[0;1];
    end
else
    if y == 1-x
        w=[.5;.5];
    elseif y < 1-x
        w=[1;0];
    else
        w=[0;1];
    end
end