function viz = neighbors_4quad(internos,coord,npt)
% routine that calculates the 05 neighboring points through MDFG - 4 quadrants
% Input data:
% npt = total number of points (internal and contour)
% coordin = array of coordinates of type [node x and z]
% internal = vector of internal points

% Operation (choice of neighbors)
% 1. First point is closest to pivot
% 2. 2nd-4th points are chosen in remaining quadrants obeying
% angle defined by tol1 tolerance
% 3. The 5th point is obtained by restarting with points discarded in item 2 and
% obeying tol2 tolerance


% loop along internal points
for w=1:size(internos,2)
    
    interno=internos(1,w);    % rated internal point
    var=find(internos==interno); % retrieves the unknown corresponding to the internal point
    
    % pivot coordinates
    xo=coord(interno,2);
    yo=coord(interno,3);

    % nearest neighbors 
    for i=1:npt
        dist(i,1)=i;
        dist(i,2)=sqrt((coord(i,2)-xo)^2+(coord(i,3)-yo)^2);
    end

    % distance vector organized
    dist_ord=sortrows(dist,2); 
    % vector of selected points for neighbors
    selec=[];
    selec(1)=dist_ord(2,1); % first point selected is closest to the pivot
    
    
    % computes the corresponding quadrant for the 'n' nearest points
    for i=2:npt
      
      ponto=dist_ord(i,1);
      
      distxr=coord(ponto,2)-xo;
      distyr=coord(ponto,3)-yo;
      
      %  Quad 1
      if distxr>0 && distyr>=0
          dist_ord(i,3)=1;
          
        % Quad 4
        else if distxr>=0 && distyr<0
              dist_ord(i,3)=4;
               
        %  Quad 2
        else if distxr<=0 && distyr>0
             dist_ord(i,3)=2;
            
        %  Quad 3
        else if distxr<0 && distyr<=0
            dist_ord(i,3)=3;
           
            end
            end
            end
      end
      
    end
      
      % definition of the original quadrants
      quad=[1 2 3 4];
      % quadrant of the nearest point selected
      quad_ponto=dist_ord(2,3);
      % removes quadrant from the nearest point of the original vector
      quad=setdiff(quad,quad_ponto);
      % sets counter for neighbors from the 2nd
      contador=2;
   
      % LOOP 1: over the neighbors to choose an additional 3 points
      % (total = 4)
      for i=3:npt
          ponto=dist_ord(i,1);      % actual point
          quad_ponto=dist_ord(i,3); % corresponding quadrant
          
          if isempty(quad)
              break
          end
                  
          % the demand for the first 4 points follows the criterion 4 quadrants
          %  includes criterion to avoid aligned points
          if ismember(quad_ponto,quad) 
            xp=coord(ponto,2);   
            yp=coord(ponto,3);   
            
            % verifies angle between straight lines formed by pivot and existing points           
            ang=[];
            for j=1:size(selec,2)
                no_ref=selec(j);          % existing selection node
                u = [xo yo 0]-[xp yp 0];  % from the pivot to the study node
                v = [xo yo 0]-[coord(no_ref,2) coord(no_ref,3) 0]; % straight from the pivot to the existing node
                ang(j) = atan2d(norm(cross(u,v)),dot(u,v));  % angle between straight lines
            end
            tol1=20; 
              if min(ang)>tol1
              selec(contador)=ponto;
              contador=contador+1;
              quad=setdiff(quad,quad_ponto);
              end
          end         
      end % end LOOP 1
      
    
      if size(selec,2)~=4
          disp('>>>>>> LOOP 1 - IMPOSSIBLE TO FIND 04 NEIGHBORHOODS <<<<<')
          disp('>>>>>>        analyze tolerance tol1           <<<<<')
          disp(['>>>>>> Point =' num2str(w)])
          
      end
      
      % removes the 4 neighbors already selected and restarts neighbors to capture the last point (5th point)
      dist_aux=dist_ord(:,1);
      dist_aux=setdiff(dist_aux,selec,'stable');
      
      % LOOP 2: choice of the last point (5th point)
      for i=1:size(dist_aux,1)
           ponto=dist_aux(i,1);     % actual point
           xp=coord(ponto,2);       
           yp=coord(ponto,3);      
           ang=[];
           for j=1:size(selec,2)
                no_ref=selec(j);          % existing selection node
                u = [xo yo 0]-[xp yp 0];  % from the pivot to the study node
                v = [xo yo 0]-[coord(no_ref,2) coord(no_ref,3) 0]; % straight from the pivot to the existing node
                ang(j) = atan2d(norm(cross(u,v)),dot(u,v));  % angle between straight lines 
           end
           tol2=20; 
              if min(ang)>tol2
              selec(5)=ponto;
              break
              end
      end  % end LOOP 2         
           
      % issues warning in the absence of 05 neighbors
      if size(selec,2)~=5
         disp('>>>>>> LOOP 2 - IMPOSSIBLE TO LOCATE 05 NEIGHBORHOODS <<<<<')
         disp('>>>>>>         analyze tolerance tol2           <<<<<')
         disp(['>>>>>> Point =' num2str(w)])
         break
      end
     
    
      % stores the selected 05 neighbors
      viz(w,:)=selec;


end    % end of the loop along the inner points

end    % end function

