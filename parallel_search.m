function[viz]=parallel_search(coord1,n,ni)
% Neighbor search by parallel programming with four cores
% input function: coord1 - coordinate matrix
% n - number of total points
% ni - number of internal points

%matlabpool open 4

viz={1 2 3 4};
parfor i=1:4
    if i==1
    viz{i}=neighbors_4quad((1:floor(ni/4)),coord1,n);
    elseif i==2
    viz{i}=neighbors_4quad((floor(ni/4)+1:floor(ni/2)),coord1,n);
    elseif i==3
    viz{i}=neighbors_4quad((floor(ni/2)+1:floor(2*ni/3)),coord1,n);
    else
    viz{i}=neighbors_4quad((floor(2*ni/3)+1:ni),coord1,n);  
    end
       
end

end