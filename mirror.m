function[x,y,nv]=mirror(x,y,ni,N,del)
% function intended to project contour points in the normal direction
% function input: x, y - coordinates of the total points
% ni - number of internal points
% N - Tilt angle tangent
% del - factor amplifier or reducer of the distance from the contour points to the virtual points
n=length(x);

aux=1;
for i=(ni+1):n
  
    xv(aux)=x(i)+N(aux,1)/del;
    yv(aux)=y(i)+N(aux,2)/del;
    aux=1+aux;
end
nv=length(xv);

x = [x xv];
y = [y yv];
