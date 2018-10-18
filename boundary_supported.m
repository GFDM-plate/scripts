function[M1]=boundary_supported(cont,i,pivo,x,y,n,ni,nf,M1,teta,viz2,D,v)

iddist=viz2(cont,1:5);

h1=x(iddist(1))-pivo(1);
h2=x(iddist(2))-pivo(1);
h3=x(iddist(3))-pivo(1);
h4=x(iddist(4))-pivo(1);
h5=x(iddist(5))-pivo(1);

k1=y(iddist(1))-pivo(2);
k2=y(iddist(2))-pivo(2);
k3=y(iddist(3))-pivo(2);
k4=y(iddist(4))-pivo(2);
k5=y(iddist(5))-pivo(2);

A=[h1 k1 h1^2/2 k1^2/2 h1*k1;...
   h2 k2 h2^2/2 k2^2/2 h2*k2;...
   h3 k3 h3^2/2 k3^2/2 h3*k3;...
   h4 k4 h4^2/2 k4^2/2 h4*k4;...
   h5 k5 h5^2/2 k5^2/2 h5*k5];

B=inv(A);

% @2/@x2 (second derivate in x) - 3a line of [B]
B_3=B(3,:);
% @2/@x2 (second derivate in y) - 4a line of [B]
B_4=B(4,:);
% (crusade derivade) - 5a line of [B]
B_5=B(5,:);

% application : d2w/dx + d2w/dy = -M/D 
[iddist_internos,ia]=setdiff(iddist,ni+1:nf); %calculate internal neighbors
M1(i+ni,n+iddist_internos)=B_3(ia)+B_4(ia); % introduction in the coefficient matrix
[iddist_virtuais,ia]=setdiff(iddist,1:n); % calculate virtual neighbors
M1(i+ni,ni+iddist_virtuais)=B_3(ia)+B_4(ia); % introduction in the coefficient matrix
M1(i+ni,i)=1/D; % generate right side of equation
   
% boundary supported: (d²w/dx² + v*d²w/dy²)cos@² + (d²w/dy² + v*d²w/dx²)sin@² + 2(1-v)d²w/(dxdy)sin@cos@ = 0
[iddist_internos,ia]=setdiff(iddist,ni+1:nf); % calculate internal neighbors    
M1(i+n,n+iddist_internos)=(B_3(ia)+v*B_4(ia))*(cos(teta(i-ni)))^2+(B_4(ia)+v*B_3(ia))*(sin(teta(i-ni)))^2+2*(1-v)*B_5(ia)*cos(teta(i-ni))*sin(teta(i-ni)); %  introduction in the coefficient matrix
[iddist_virtuais,ia]=setdiff(iddist,1:n); % calculate virtual neighbors
M1(i+n,ni+iddist_virtuais)=(B_3(ia)+v*B_4(ia))*(cos(teta(i-ni)))^2+(B_4(ia)+v*B_3(ia))*(sin(teta(i-ni)))^2+2*(1-v)*B_5(ia)*cos(teta(i-ni))*sin(teta(i-ni)); % introduction in the coefficient matrix


end