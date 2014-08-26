function [bdofQ2,nbnodes]=ComputeBoundaryDOFQuad(mapping,VBasis,ZSTRIDE)

NZ=size(mapping,2);
left=[];
right=[];
top=[];
bot=[];
switch func2str(VBasis)
    case {'Q1Basis', 'Q1bBasis'}
        for i=1:ZSTRIDE:NZ
            left=[left,mapping(1,i),mapping(4,i)];
        end

        for i=ZSTRIDE:ZSTRIDE:NZ
            right=[right,mapping(2,i),mapping(3,i)];
        end

        for i=NZ-ZSTRIDE+1:NZ
            top=[top,mapping(4,i),mapping(3,i)];
        end

        for i=1:ZSTRIDE
            bot=[bot,mapping(1,i),mapping(2,i)];
        end
    case 'Q2Basis'
        for i=1:ZSTRIDE:NZ
            left=[left,mapping(1,i),mapping(8,i),mapping(4,i)];
        end

        for i=ZSTRIDE:ZSTRIDE:NZ
            right=[right,mapping(2,i),mapping(6,i),mapping(3,i)];
        end

        for i=NZ-ZSTRIDE+1:NZ
            top=[top,mapping(4,i),mapping(7,i),mapping(3,i)];
        end

        for i=1:ZSTRIDE
            bot=[bot,mapping(1,i),mapping(5,i),mapping(2,i)];
        end
end

bdofQ2=[bot,right,top,left]';
nbnodes=[length(bot),length(right),length(top),length(left)];