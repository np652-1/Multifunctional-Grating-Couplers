function[dist] = getDistVertLineSegment(xy,xLine,yLine)
% Get the distance between a point and a vertical line
if (xy(2) > yLine(2)) % above top of wall
    dist = norm(xy - [xLine(2) yLine(2)]);
elseif (xy(2) < yLine(1)) % below bottom of wall
    dist = norm(xy - [xLine(1) yLine(1)]);
else % within the y-bounds of the wall
    dist = abs(xLine(2) - xy(1));
end
end