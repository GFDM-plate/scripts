function[teta,N]=normal_contour(coord,conec)
% Function input: cood - total domain coordinates
% connec - contour connectivity
% function output: theta - inclination angle from normal to contour with respect to axis
% horizontal
% N - Normal direction at the contour.

Vertices=[coord(:,2) coord(:,3)];
Lines=[conec(:,2) conec(:,3)];

N=LineNormals2D(Vertices,Lines);

N1=N(:,1);
N2=N(:,2);

teta = atan(N(:,2)./N(:,1));
teta(isnan(teta)) = [];

N1(isnan(N1)) = [];
N2(isnan(N2)) = [];
N=[N1 N2];
end