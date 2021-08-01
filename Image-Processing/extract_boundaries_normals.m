function [ptCloud_IN, ptCloud_OUT] = extract_boundaries_normals(points, n_theta_points, min_wall_thickness)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script takes in a point cloud of the artery and delineates the inner
% and outer boundaries along with their normals.
% INPUT :   points            (Nx3) - 3D points 
%           n_theta_points      (1x1) - number of theta points to be sampled in each slice
%           min_wall_thickness  (1x1) - minimum wall thickness in case the outer wall needs to be pushed away
% OUTPUT:   ptCloud_IN          (pointCloud object) - pointcloud of the inner wall WITH normals           
%           ptCloud_OUT         (pointCloud object) - pointcloud of the outer wall WITH normals     
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin <2)    n_theta_points = 100;  end % points per slice in inner and outer boundaries
if(nargin <3)    min_wall_thickness = 0; end 

points = sortrows(points,3);        % Sort the coordinates by the z value
unique_z = unique(points(:,3));     % Get all unique z values (all the slices)
OUT = []; IN = []; normals_OUT = []; normals_IN = []; 

% -- LOOP THROUGH ALL THE SLICES -- %
for iSlice = 1:length(unique_z)
    if(mod(unique_z(iSlice),0.4)) > 0
        continue;
    end    
    % Get current slice coordinates only
    slice = points(points(:,3) == unique_z(iSlice),:);
    [theta,rho] = cart2pol(slice(:,1),slice(:,2));            % Convert to polar coordinates
    
    % -- FOR EACH DEGREE INTERVAL, FIND THE INNERMOST AND OUTERMOST POINTS -- %
    slice_polar = [];  slice_polar(:,1) = theta;  slice_polar(:,2) = rho;  slice_polar(:,3) = 1:length(theta);
    slice_polar = sortrows(slice_polar,1);
    interval = pi/n_theta_points;   
    start = -pi; stop = start + interval;
    i = 1;
    
    while stop <= pi
        
        mini_slice = slice_polar(slice_polar(:,1)>=start & slice_polar(:,1)<=stop,:);
        
        [~,j] = max(mini_slice(:,2)); % Find point with the highest radius in the current interval
        outer_bound(i) = mini_slice(j,3);
        [~,j] = min(mini_slice(:,2)); % Find the points with the lowest radius in the current interval
        inner_bound(i) = mini_slice(j,3);
        
        start = start + interval; stop = stop + interval;        i = i+1 ;
        
    end
    
    OUT = [OUT;slice(outer_bound,:)];
    IN = [IN;slice(inner_bound,:)];
    
    slice_outer = slice(outer_bound,:); slice_inner = slice(inner_bound,:);
    normals_OUT = [normals_OUT;get_normals(slice_outer)];
    normals_IN = [normals_IN;get_normals(slice_inner)];
    
end

% Push the outer points if they are too close to the inner points
if(min_wall_thickness)
    OUT = push_boundaries(IN,OUT,min_wall_thickness);
end

% ---------- CREATE POINT CLOUDS AND COMPUTE NORMALS ---------------
ptCloud_OUT = pointCloud([OUT]);
ptCloud_OUT.Normal = normals_OUT;

ptCloud_IN = pointCloud([IN]);
ptCloud_IN.Normal = normals_IN;


end


