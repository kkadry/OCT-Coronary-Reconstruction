%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script takes in a point cloud of the artery and delineates the inner
% and outer boundaries along with their normals, as well as each in-frame plaque component.
% INPUT :   plaque           (Nx4) txt file - 3D points and material assignment number for multiple frames 
%         
% OUTPUT:   Frame          (matlab dataframe) - 3D boundary point data for vessel borders and each plaque component          
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
addpath(genpath(pwd))

%% Import point cloud data
dataname="";
plaque=importdata(dataname);
points=plaque(:,1:3);
disp('Plaque Data Imported')

%% Extract inner and outer boundaries from point cloud
[PCin,PCout]=extract_boundaries_normals(points,100,0.5);

%% Extract discrete points of inner and outer boundaries
pcin=PCin.Location;
pcout=PCout.Location;
normalsin=PCin.Normal;
normalsout=PCout.Normal;
disp('Boundary Point Clouds Obtained')


%% Loop over all frames
dx = 0.022; %OCT intra-frame resolution
z = unique(pcout(:,3));
z = unique(points(:,3));
dz = z(2)-z(1); %OCT inter-frame resolution, mm
frames = 1:length(z);
slices = [];
luminal_lipid_frame=[];

for frame = frames


    % Identify plaque points in frame of interest
    idx = find(points(:,3) == z(frame));
    for i = 1:5
        plaque_frame(:,i) = plaque(idx,i);
    end
    points_frame = plaque_frame(:,1:3);
   
    
    sz = max(size(unique(points(:,1)),1),size(unique(points(:,2)),1))+20; %pad with 20
    x = round((plaque_frame(:,1)-min(points(:,1)))./dx+10); %add 10 to avoid border
    y = round((plaque_frame(:,2)-min(points(:,2)))./dx+10); %add 10 to avoid border
    
    slice = zeros(sz, sz);
    sliceID = slice;
   
    for n = 1:size(plaque_frame,1)
        slice(x(n),y(n)) = plaque_frame(n,4); %slice is the annotated image in selected frame
    end
    slice_orig=slice;
    %Filtering small parts (except for calcium and lipid)    
    thresh = 300*[0,1,0,1,1];
    for i = 1:4
        bin_slice = double(slice == i);
        bin_slice_thresh = double(bwareaopen(bin_slice,thresh(i)));
        removed = bin_slice - bin_slice_thresh;
        if sum(removed(:)) > 0
            removed_comp = bwconncomp(removed);
            for n = 1:removed_comp.NumObjects
                removed_n = zeros(size(removed));
                removed_n(removed_comp.PixelIdxList{n}) = 1;
                dilatedImage = imdilate(removed_n, true(3));
                borderPixels = bwperim(dilatedImage);
                neighbor = borderPixels .* slice;
                v = unique(neighbor(neighbor > 0));

                if size(v,1) == 1
                    slice(removed_n > 0) = v;

                else
                    s = [];
                    sID=[];
                    for m = 1:size(v)
                        s(m) = sum(neighbor(:) == v(m));
                    end
                    [~,pos] = max(s);
                    slice(removed_n > 0) = v(pos);

                end
            end
        end
    end
    

%Redoing the connectivity analysis
sliceID=zeros(size(sliceID));
i2=unique(slice)';
for i=i2(2:end)
    mask=double(slice==i);
    bw=bwconncomp(mask);
    for j=[1:bw.NumObjects]
        sliceID(bw.PixelIdxList{j})=j;        
    end 
end


%Adding a minimum thickness for the luminal lipid
growfac=5;
lumask=double(slice==0);

%Selecting the center by getting the size of the frame
framecenter=size(slice)/2;
%Getting the innerlum connected component
innerlum=bwselect(lumask,framecenter(1),framecenter(2));
%Setting innerlumdilate as innerlum so i can dilate it freely
innerlumdilate=innerlum;
%Dilating the inner lumen component to find intersection with lipid
innerlumdilate=imdilate(innerlumdilate, true(3));

%Take whatever plaque component intersects with dilated lipid
bordertissue=(innerlumdilate-innerlum).*slice;
%Take only lipid from that
lumenlipid=double(bordertissue==3);
lumenlipidilate=lumenlipid;
%Dilate the luminal lipid and keep taking the intersection with normal
%lipid and dilate that again
for i=[1:growfac]
%Dilating the lumlip
lumenlipidilate=imdilate(lumenlipidilate, true(3));
%Taking the intersection with lipid
lip=double(slice==3);
lumenlipidilate=lip.*lumenlipidilate;
end

%Combine it with innerlum
lumpluslip=(lumenlipid+innerlum);
%conmask=bwselect(lumpluslip,160,160)-innerlum;
%conmask=bwselect(lumpluslip,framecenter(1),framecenter(1))-innerlum;
conmask=lumenlipidilate;
    if sum(unique(conmask))>0
    close all
    lipidx=find(conmask==1);
    slice(lipidx)=4;
    %Updating the luminal lipid vector to keep track
    luminal_lipid_frame(end+1)=frame;
    end


%Getting the distance transform to find thickness of fibrous cap
[Dmap,Didx]=bwdist(innerlum);
%Finding the lumen centroid and area
lumprop=regionprops(innerlum,'Area','Centroid');
%Adjusting the units
%Multiplying the area by area of a pixel
lumprop.Area=lumprop.Area*dx^2;
%Adjusting the centroid from pixel to real space
lumprop.Centroid(1)=(lumprop.Centroid(1)-158.5)*dx;
lumprop.Centroid(2)=(lumprop.Centroid(2)-158.5)*dx;
    %% Identify material borders in slice
    
    lF = []; %initialize faces
    V = []; %initialize vertices
    
    % Initialize Frame-structure
    slice_labels = unique(nonzeros(slice));
    if ~ismember(1,slice_labels)
        Frame{frame}.Mat1 = [];
    end
    if ~ismember(2,slice_labels)
        Frame{frame}.Mat2 = [];
    end
    if ~ismember(3,slice_labels)
        Frame{frame}.Mat3 = [];
    end
    if ~ismember(4,slice_labels)
        Frame{frame}.Mat4 = [];
    end
    if ~ismember(5,slice_labels)
        Frame{frame}.Mat5 = [];
    end
    
    iter1 = 0; iter2 = 0; iter3 = 0; iter4 = 0; iter5 = 0;
    for mat = [1,2,3,4,5] %loop over all labels
        label = (slice == mat);
        labelIDs = unique(nonzeros(label .* sliceID));
        if iscolumn(labelIDs)
            labelIDs = labelIDs';
        end
        for labelID = labelIDs
            if sum(label(:)) == 0
                continue
            end
            label = ((sliceID == labelID) .* (slice == mat));
            if mat~=3
                % Dilate outwards
                label = imdilate(label, ones(3,3));
                %If it's fibrous (mat=4) then dilate twice again to make intersections
                %easier for boolean operations
                if mat==4
                    label=imdilate(label, ones(3,3));
                    label=imdilate(label, ones(3,3));
                end
            end
            
            
            %Storing a flat version of label for property extraction
            labelflat=label;
            label(:,:,1) = label;
            label(:,:,2) = label;
            
            if ~exist('X')
                [X,Y,Z] = meshgrid(1:size(label,2),1:size(label,1),1:size(label,3));
                cutoff = 0.5; %define cut-off for isosurface extraction
            end

            % Isosurface - determine surface along boundary
            [f,v] = isosurface(X,Y,Z,label,cutoff);

            %Collapse 3rd dimension
            v(:,3) = 1;

            %Remove redundant/duplicate patches
            [v, f] = patchslim(v, f);

            % Convert back from pixel to real dimesions
            v(:,1) = (v(:,1)-10).*dx+min(points(:,1));
            v(:,2) = (v(:,2)-10).*dx+min(points(:,2));
            v(:,3) = z(frame);
            
            %smooth boundaries
            FV.vertices = v;
            FV.faces = f;
            FV = smoothpatch(FV,1,5,1);

            %split separate boundaries
            fvOut = splitFV(FV.faces, FV.vertices);
            for j = 1:size(fvOut,1)
                FV.vertices = fvOut(j).vertices;
                FV.faces = fvOut(j).faces;
                
                % sort data
                [~,pos] = min(FV.vertices(:,1));
                v_new = FV.vertices(pos,:);
                f_new = 1;
                FV.vertices(pos,:) = [];
                for i = 1:size(FV.vertices,1)
                    pos = dsearchn(FV.vertices,v_new(i,:));
                    v_new = [v_new; FV.vertices(pos,:)];
                    f_new = [f_new; i+1];
                    FV.vertices(pos,:) = [];

                end
                % repeat last data
                v_new = [v_new; v_new(1,:)];
                f_new = [f_new; 1];
                FV.vertices = v_new;
                FV.faces = f_new;
                
                %Finding region properties
                geoprop=regionprops(labelflat,'Area','Centroid','MajorAxisLength','MinorAxisLength','Orientation');
                %Adjusting the units
                %Multiplying the centroid by pixel area
                geoprop.Area=geoprop.Area*dx^2;
                %Adjusting centroid position to real space
                geoprop.Centroid(1)=(geoprop.Centroid(1)-158.5)*dx;
                geoprop.Centroid(2)=(geoprop.Centroid(2)-158.5)*dx;
                %Multiplying the axis lengths with pixel length
                geoprop.MajorAxisLength=geoprop.MajorAxisLength*dx;
                geoprop.MinorAxisLength=geoprop.MinorAxisLength*dx;
                %Orientation is unchanged
                
                %Finding the distance to lumen
                %Distance transform
                labelDmap=Dmap.*labelflat;
                %Taking out the 0 distance elements
                labelDmap(labelDmap==0)=100;
                %Finding the minimum distance to lumen
                minthick=double(min(min(labelDmap)));
                %Storing the value while correcting the units (in um)
                geoprop.Min_thickness=minthick*dx*1000;
                % Append lists
                if mat == 1
                    iter1 = iter1+1;
                    Frame{frame}.Mat1{iter1}.vertices = FV.vertices;
                    Frame{frame}.Mat1{iter1}.faces = FV.faces;
                    Frame{frame}.Mat1{iter1}.ID = labelID;
                    Frame{frame}.Mat1{iter1}.geom=geoprop;
                elseif mat == 2
                    iter2 = iter2+1;
                    Frame{frame}.Mat2{iter2}.vertices = FV.vertices;
                    Frame{frame}.Mat2{iter2}.faces = FV.faces;
                    Frame{frame}.Mat2{iter2}.ID = labelID;
                    Frame{frame}.Mat2{iter2}.geom = geoprop;
                elseif mat == 3
                    iter3 = iter3+1;
                    Frame{frame}.Mat3{iter3}.vertices = FV.vertices;
                    Frame{frame}.Mat3{iter3}.faces = FV.faces;
                    Frame{frame}.Mat3{iter3}.ID = labelID;
                    Frame{frame}.Mat3{iter3}.geom = geoprop;

                elseif mat == 4
                    iter4 = iter4+1;
                    Frame{frame}.Mat4{iter4}.vertices = FV.vertices;
                    Frame{frame}.Mat4{iter4}.faces = FV.faces;
                    Frame{frame}.Mat4{iter4}.ID = labelID;
                    Frame{frame}.Mat4{iter4}.geom = geoprop;
                else 
                    iter5 = iter5+1;
                    Frame{frame}.Mat5{iter5}.vertices = FV.vertices;
                    Frame{frame}.Mat5{iter5}.faces = FV.faces;
                    Frame{frame}.Mat5{iter5}.ID = labelID;
                end
            end
            clear FV fvOut
        end
    end    
    Frame{frame}.slice = slice;
       
    % Identify inner/outer borders in each slice
    label = (slice ~= 0);
   
    % Identify inner/outer borders in each slice
    label = (slice ~= 0);
    label(:,:,1) = label;
    label(:,:,2) = label;
    [X,Y,Z] = meshgrid(1:size(label,2),1:size(label,1),1:size(label,3));
    cutoff = 0.5; %define cut-off for isosurface extraction
    
    % Isosurface - determine surface along boundary
    [f,v] = isosurface(X,Y,Z,label,cutoff);

    %Collapse 3rd dimension
    v(:,3) = 1;

    %Remove redundant/duplicate patches
    [v, f] = patchslim(v, f);

    % Convert back from pixel to real dimesions
    v(:,1) = (v(:,1)-10).*dx+min(points(:,1));
    v(:,2) = (v(:,2)-10).*dx+min(points(:,2));
    v(:,3) = z(frame);
            
    %smooth boundaries
    FV.vertices = v;
    FV.faces = f;
    FV = smoothpatch(FV,1,5,0.75);
    

    fvOut = splitFV(FV.faces, FV.vertices);
    for j = 1:size(fvOut,1)
        FV.vertices = fvOut(j).vertices;
        FV.faces = fvOut(j).faces;

        % sort data
        [~,pos] = min(FV.vertices(:,1));
        v_new = FV.vertices(pos,:);
        f_new = 1;
        FV.vertices(pos,:) = [];
        for i = 1:size(FV.vertices,1)
            pos = dsearchn(FV.vertices,v_new(i,:));
            v_new = [v_new; FV.vertices(pos,:)];
            f_new = [f_new; i+1];
            FV.vertices(pos,:) = [];
        end
        % repeat last data
        v_new = [v_new; v_new(1,:)];
        f_new = [f_new; 1];
        if j == 1
            FV1.vertices = v_new;
            FV1.faces = f_new;
        elseif j == 2
            FV2.vertices = v_new;
            FV2.faces = f_new;
        end       
    end


    %Determine which is inner/outer
    A1 = polyarea(FV1.vertices(:,1),FV1.vertices(:,2));
    A2 = polyarea(FV2.vertices(:,1),FV2.vertices(:,2));
    if A1 < A2
        %Setting the minimum wall thickness
        min_wall_thickness=0.5;
        %Push out the boundaries of the outer lumen
        Frame{frame}.Inner.vertices = FV1.vertices;
        Frame{frame}.Inner.faces = FV1.faces;
        Frame{frame}.Outer.vertices = FV2.vertices;
        Frame{frame}.Outer.faces = FV2.faces;
        
        Frame{frame}.Inner.geom =lumprop;
    else
        Frame{frame}.Inner.vertices = FV2.vertices;
        Frame{frame}.Inner.faces = FV2.faces;
        Frame{frame}.Outer.vertices = FV1.vertices;
        Frame{frame}.Outer.faces = FV1.faces;
        
        Frame{frame}.Inner.geom =lumprop;
    end
    
    
    fprintf('Slice %i complete out of %i\n',frame,length(frames));
    clear plaque_frame X Y Z
end
%Saving all frame indices where luminal lipid was found
pname="PatientMat.mat";
save(pname,'Frame');
%% OUTPUT EXPLANATION %%
%
% Output is given in a cell array form, with the size given by the number of 
% frames in the OCT pullback
% 
% in frame_i:
% 
% Frame{i}.Mat1 - structure list of cells, each entity describing a closed boundary, e.g.
%     Frame{i}.Mat1{1}:
%         .vertices    - vertice points (nodal points)
%         .faces       - connectivity list
%         .ID          - ID of component (taken from Max's data)
% Frame{i}.Mat2 - same as above
% Frame{i}.Mat3 - same as above
% Frame{i}.Mat4 - same as above
% Frame{i}.Mat5 - same as above
% Frame{i}.slice - matrix of annotated label, given in pixel idx
% Frame{i}.Inner - structure for inner detection, containing
%     .vertices   - vertice points (nodal points)
%     .faces      - connectivity list
% Frame{i}.Outer - structure for outer detection, containing
%     .vertices   - vertice points (nodal points)
%     .faces      - connectivity list
%
