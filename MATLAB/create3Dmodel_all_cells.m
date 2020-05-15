% This is a version of the script we used to convert the ExSeq transcripts
% and morphology traced in VAST-Lite into objects that can be loaded into
% the ExSeqViewer.
% This is research-grade code that has worked for an N=1 computers so far,
% so please reach out to dgoodwin@mit.edu if you would like to learn more
% about getting this running for your own experiments! Best, Dan

addpath('3D graphics/'); 

%In the ExSeq paper, we traced three fields of view which were numbered
%8-10
fovs = 8:10;


loadParameters;
RESULTSDIR = ''; %Where is the transcript_object file?
OFDIR = ''; %Where is the /bin/data folder of your openframeworks app?/

if ~exist('transcript_objects','var')
   load(fullfile(RESULTSDIR,'transcript_objects.mat')); 
end

%Taking only a fraction of the morphology points to create the mesh: this
%helps keep the file sizes reaosnable, and also to smooth out the
%morphology
DOWNSAMPLE_PERCENTAGE = 0.8;

%When we traced neurons in VAST-lite, each part of the morphology 
%(soma, nucleus and dendrite) were given unique IDs. We then numbered each
%cell, and manually noted the cell number - soma id and dendrite below in
%this data structure:

cellparts{8} = [3,6;
    2,1;
    7,4;
    12,9;
    8,5;];
cellparts{9} = [ 13, 8;
    16, 9;
    1, 2;];
cellparts{10} = [8,15;
    4,13;
    16,7;
    2,1;];

%For each field of view, load the segmented file
for f_idx = 1:length(fovs)
    
    fov_index = fovs(f_idx);
    
    segmentation = load3DTif_uint16(fullfile(RESULTSDIR,sprintf('xy%.2i_segmentation.tif',fov_index)));
    
    % Downsample it for speed
    seg_downsample = imresize3(segmentation,1/params.DOWNSAMPLE_RATE,'nearest');
    

    %The segmentation are complete volumes, but we just want the surface
    %for the 3D visualizer. We use the gradient to get the surface.
    %If we want to piece together complete cells, use these two lines
    for cell_idx = 1:size(cellparts{fov_index},1)
        [Gx,Gy,Gz] = imgradientxyz((seg_downsample==cellparts{fov_index}(cell_idx,1)) | (seg_downsample==cellparts{fov_index}(cell_idx,2)));
        
        mag = Gx.^2+Gy.^2+Gz.^2;
        [e2, thresh2] = canny(mag);
        
        
        
        point_indices = find(e2);
        [x,y,z] = ind2sub(size(seg_downsample),point_indices);
        
        %downsample randomly:
        keepers = randperm(length(x),round(DOWNSAMPLE_PERCENTAGE*length(point_indices)));
        
        x = x(keepers);
        y = y(keepers);
        z = z(keepers);
        
        %Hardcoding the 2048x204x256 image volume size into a 2x downsample
        %coordinates, getting a center poitn
        morphmean_x = 1024/2; 
        morphmean_y = 1024*f_idx-1024/2; 
        morphmean_z = 128/2;

        %remove the mean:
        x = x - morphmean_x;
        y = y - morphmean_y;
        z = z - morphmean_z;
        p = [x y z];
        if numel(p)==0 %In the edge case where the morphology index has nothing
            continue
        end
        % Run  program to build a shape
        [t]=MyCrustOpen(p);
        
        plyWrite(p,t,fullfile(OFDIR,sprintf('fov_%i_neuron_%i.ply',fov_index, cell_idx)));
        
        fprintf('Completed morphology %i\n',cell_idx);
    end
    
    fprintf('Done with morphology for FoV %i!\n',fov_index);
    
    
    % Now go through the transcript objects file to save all reads in the
    % cells into a CSV for the OpenFrameworks app to load
    
    fov_indices = cell2mat(cellfun(@(x) [x.fov ==fov_index]',transcript_objects,'UniformOutput',0));
    
    transcript_objects_active = transcript_objects(fov_indices);
    
    punctapos = cell2mat(cellfun(@(x) [x.pos]',transcript_objects_active,'UniformOutput',0));
    punctapos = punctapos';
    
    punctapos_orig = punctapos;
    
    %the matlab centroid switches the x and y,
    punctapos(:,1) = punctapos_orig(:,2)/params.DOWNSAMPLE_RATE-morphmean_x;
    punctapos(:,2) = punctapos_orig(:,1)/params.DOWNSAMPLE_RATE-morphmean_y;
    punctapos(:,3) = punctapos_orig(:,3)/params.DOWNSAMPLE_RATE-morphmean_z;
    
    %Use a numbering system to map the readtypes to an integer to simplify the
    %C++ side in OpenFrameworks
    readtypes = {'intron','exon','intergenic','rrna','noalignment','refseqnotingenome'};
    
    output_file = fullfile(OFDIR,sprintf('fov_%i_puncta.csv',fov_index));
    
    fid = fopen(output_file,'w');
    for idx= 1:size(punctapos,1)
        tobj = transcript_objects_active{idx};

        gene_symbol = '';
        if isfield(tobj,'genename')
            gene_symbol = tobj.genename;
        else
            gene_symbol = tobj.readtype;
        end
        gene_symbol = replace(gene_symbol,',','');
        if strcmp(gene_symbol,'ambiguous')
            continue;
        end
        tobj_readtype = tobj.readtype;
        
        %Make a manual correction to the visualizer: about 25% of intergenics
        %actually had a gene symbol and aligned to refseq, meaning there was
        %was something funky with the refseq gene's alignment to the genome
        if (strcmp('intergenic',tobj_readtype) && isfield(tobj,'genename'))
            tobj_readtype = 'exon';
        end
        readtype = find(strcmp(readtypes,tobj_readtype)); %returns a cell array of hits
        didalign = 1;
        incell = isfield(tobj,'cell_number');
        
        
        fprintf(fid,'%i,%i,%i,%s,%i,%i,%i\n', round(punctapos(idx,1)),round(punctapos(idx,2)),round(punctapos(idx,3)),...
            gene_symbol,readtype,didalign,incell );
        
    end
    fprintf('For loop complete\n');
    fclose(fid);
    fprintf('Done!\n')
end

