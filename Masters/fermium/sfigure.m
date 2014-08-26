function h = sfigure(varargin)
% SFIGURE  Create figure window (minus annoying focus-theft).
%
% Usage is identical to figure.
%
% Daniel Eaton, 2005
%
% See also figure

if nargin>=1 
    h=varargin{1};
	if ishandle(h)
		set(0, 'CurrentFigure', h);
	else
		h = figure(h);
    end
else
	h = figure;
end

if nargin>=2
    set(h,varargin{2:end});
end