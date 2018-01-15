close all; clear all; clc;
xL = 0; yL = 0; xM = 0.5; % defin coords
yM = 1; xN = 2; yN = 0; % define coords
AL = 10; AM = 15; AN = 1; % define function

xv = [xL xM xN xL]; 
yv = [yL yM yN xL]; 
Av = [AL AM AN AL]; % for plotting

figure
plot3(xv, yv,  Av,'-o'); 
grid on; hold on
xlabel('x'); 
ylabel('y'); 
zlabel('Function Value: A')

Mxy = [ 1 xL yL; 1 xM yM; 1 xN yN];
Mab = inv(Mxy);  % equation 4.9
aL = Mab(1,1); aM = Mab(1,2); aN = Mab(1,3);
bL = Mab(2,1); bM = Mab(2,2); bN = Mab(2,3);
cL = Mab(3,1); cM = Mab(3,2); cN = Mab(3,3);

xPlot = linspace(min([xL xM xN]), max([xL xM xN]),80);
yPlot = linspace(min([yL yM yN]), max([yL yM yN]),90);

for xIndex = 1:length(xPlot)
    for yIndex = 1:length(yPlot)
        x = xPlot(xIndex);
        y = yPlot(yIndex);
        if inpolygon(x,y,xv,yv);
            A(yIndex,xIndex) = AL*(aL+bL*x+cL*y)+ ...
                AM*(aM+bM*x+cM*y)+AN*(aN+bN*x+cN*y);
        else
            A(yIndex,xIndex) = nan;
        end
    end
end
mesh(xPlot, yPlot, A)



