% Stores an active plot as a 3D model for printing

a = get(gca,'Children');
xdata = get(a, 'XData');
ydata = get(a, 'YData');
zdata = get(a, 'ZData');

s = surf(xdata, ydata, zdata ,'FaceAlpha',0.2 );
s.EdgeColor = 'none';

stlwrite('test.stl',xdata,ydata,zdata,'mode','ascii')