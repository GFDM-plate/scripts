% FEDERAL UNIVERSITY OF PERNAMBUCO - UFPE
% PROGRAM USING MDFG FOR SOLUTION OF PLATE FLEXIBLE PROBLEM
% Case 04 – irregular plate with mixed boundary (clamped and simply supported)
% AUTHOR: Ferreira., A. C. A

clear
clc
close all
format SHORTG

%% 1. Preparation of the data

%----------------------------INPUT-----------------------------------------
% Material characteristics 
E=21287000;
v=0.2;

% transverse load
q=100;

% thickness 
h=0.8;

% factor amplifier or reducer of distance from contour points to virtual points
ft=.8;

% (x,y) - coordinates of points (internal + contour)
% coord - coordinates exported from Gid
% conec - connectivity in the contour do Gid
% pointsgid - numbering of points from Gid
[x,y,nc,coord,conec,pointsgid]=pre_process;
%--------------------------------------------------------------------------

% flexural rigidity of the plate
D=E*h^3/(12*(1-v^2));

n=length(x);  % number of total points
ni=n-nc;  % number of internal points

M1=sparse(2*ni+nc); % creation of the coefficient matrix
coord1=[(1:n)' x' y']; % matrix with numbering of coordinates
 
[viz]=parallel_search(coord1,n,ni); % parallel search of neighbors (4 cores)
viz1=[viz{1};viz{2};viz{3};viz{4}]; % store all neighbors

%% 2. Assembly of the stiffness

% Loop of internal points (variables M and w)
for i=1:ni
       
    pivot=[x(i) y(i)]; % pivot point

    iddist=viz1(i,1:5); % 05 neighbors of the pivot i
    
    h1=x(iddist(1))-pivot(1);
    h2=x(iddist(2))-pivot(1);
    h3=x(iddist(3))-pivot(1);
    h4=x(iddist(4))-pivot(1);
    h5=x(iddist(5))-pivot(1);

    k1=y(iddist(1))-pivot(2);
    k2=y(iddist(2))-pivot(2);
    k3=y(iddist(3))-pivot(2);
    k4=y(iddist(4))-pivot(2);
    k5=y(iddist(5))-pivot(2);

    A=[h1 k1 h1^2/2 k1^2/2 h1*k1;...
   h2 k2 h2^2/2 k2^2/2 h2*k2;...
   h3 k3 h3^2/2 k3^2/2 h3*k3;...
   h4 k4 h4^2/2 k4^2/2 h4*k4;...
   h5 k5 h5^2/2 k5^2/2 h5*k5]; 

    B=inv(A); % equation 17.b

% @2/@x2 (second derivate in x) - 3a line of [B]
    B_3=B(3,:);
% @2/@y2 (second derivate in y) - 4a line of [B]
    B_4=B(4,:);

% d2M/dx + d2M/dy = -q (contour points are unknown)
    M1(i,iddist)=B_3+B_4; % introduction in the coefficient matrix
    M1(i,i)=-sum(B_3+B_4); % generate main diagonal
% d2w/dx + d2w/dy = -M/D (w=0 in contour)
   [iddist_internos,ia]=setdiff(iddist,ni+1:n); % calculate internal neighbors
   M1(i+ni,n+iddist_internos)=B_3(ia)+B_4(ia); % introduction in the coefficient matrix
   M1(i+ni,i+n)=-sum(B_3+B_4); % generate main diagonal
   M1(i+ni,i)=1/D; % generate right side of equation

end

[teta,N]=normal_contour(coord,conec); % Normal directions at the contour (see Fig.1)

[x,y,nv]=mirror(x,y,ni,N,ft); % creation of the virtual points (projection of the points in the normal direction of the contour)

nf=n+nv; % total number of points (including virtual)
coord2=[(1:nf)' x' y']; % matrix with numbering of the new coordinates
[viz2]=neighbors_4quad((ni+1:n),coord2,nf); % calculation of neighbors in the contour

% loop at the boundary points (variable w) and applying the boundary condition
cont=1;
for i=(ni+1):n  
    pivot=[x(i) y(i)]; % pivot point
    if pivot(2)<18.84   
        [M1]=boundary_clamped(cont,i,pivot,x,y,n,ni,nf,M1,teta,viz2,D); % boundary clamped    
        cont=cont+1; 
    else   
        [M1]=boundary_supported(cont,i,pivot,x,y,n,ni,nf,M1,teta,viz2,D,v); % boundary supported 
        cont=cont+1; 
    end
end

%% 3. Solution of the system of linear equation

M2=sparse(2*ni+nc+nv,1);
M2(1:ni,1)=-q; % loading matrix

sol=M1\M2; % solution
sol=full(sol); % full format
w=cat(1,sol(n+1:2*ni+nc),zeros(nc,1)); % displacements of all points of plate (without virtual points)

%% 4. Post Processing

disp(['maximum displacement [m]: ', num2str(max(w))]) 
deslocy=zeros(n,1);
sol=[pointsgid deslocy deslocy -w]; % results to be exported to the GID

% Calculation of bending moments (internal loop)
X=2*ni+nc+nv; % number of unknowns of the problem
w=cat(2,w',sol(2*ni+nc+1:X)); % displacements of all points of plate (with virtual points)
for i=1:ni
    
    pivot=[x(i) y(i)]; % pivot point
      
    iddist=viz1(i,:); % 05 neighbors of the pivot i
      
    h1=x(iddist(1))-pivot(1);
    h2=x(iddist(2))-pivot(1);
    h3=x(iddist(3))-pivot(1);
    h4=x(iddist(4))-pivot(1);
    h5=x(iddist(5))-pivot(1);
    
    k1=y(iddist(1))-pivot(2);
    k2=y(iddist(2))-pivot(2);
    k3=y(iddist(3))-pivot(2);
    k4=y(iddist(4))-pivot(2);
    k5=y(iddist(5))-pivot(2);
    
    A=[h1 k1 h1^2/2 k1^2/2 h1*k1;
   h2 k2 h2^2/2 k2^2/2 h2*k2;
   h3 k3 h3^2/2 k3^2/2 h3*k3;
   h4 k4 h4^2/2 k4^2/2 h4*k4;
   h5 k5 h5^2/2 k5^2/2 h5*k5]; 
   
   f=[w(iddist(1))-w(i);w(iddist(2))-w(i);w(iddist(3))-w(i);w(iddist(4))-w(i);w(iddist(5))-w(i)];
   df=A\f; % equation 17a
       
   Mx(i)=-D*(df(3)+v*df(4)); % bending moments in relation to the x-axis
   My(i)=-D*(df(4)+v*df(3)); % bending moments in relation to the y-axis
   Mxy(i)=D*(1-v)*df(5); % twisting moment 
  
end

% Calculation of bending moments (loop in contour)
cont=0;
for i=(ni+1):n
   
    cont=cont+1;
    pivot=[x(i) y(i)]; % pivot point
       
    iddist=viz2(cont,:); %  05 neighbors of the pivot i 
 
    h1=x(iddist(1))-pivot(1);
    h2=x(iddist(2))-pivot(1);
    h3=x(iddist(3))-pivot(1);
    h4=x(iddist(4))-pivot(1);
    h5=x(iddist(5))-pivot(1);

    k1=y(iddist(1))-pivot(2);
    k2=y(iddist(2))-pivot(2);
    k3=y(iddist(3))-pivot(2);
    k4=y(iddist(4))-pivot(2);
    k5=y(iddist(5))-pivot(2);

   A=[h1 k1 h1^2/2 k1^2/2 h1*k1;
   h2 k2 h2^2/2 k2^2/2 h2*k2;
   h3 k3 h3^2/2 k3^2/2 h3*k3;
   h4 k4 h4^2/2 k4^2/2 h4*k4;
   h5 k5 h5^2/2 k5^2/2 h5*k5];

   f=[w(iddist(1))-w(i);w(iddist(2))-w(i);w(iddist(3))-w(i);w(iddist(4))-w(i);w(iddist(5))-w(i)];
   df=A\f; % equation 17a
  
   Mx(i)=-D*(df(3)+v*df(4)); % bending moments in relation to the x-axis
   My(i)=-D*(df(4)+v*df(3)); % bending moments in relation to the y-axis
   Mxy(i)=D*(1-v)*df(5); % twisting moment 
   
end

sol_Mx=[pointsgid Mx']; % bender results to be exported to the GID
sol_My=[pointsgid My']; % bender results to be exported to the GID
sol_Mxy=[pointsgid Mxy']; % bender results to be exported to the GID
    


