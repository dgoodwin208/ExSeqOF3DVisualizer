%taken from https://www.mathworks.com/matlabcentral/fileexchange/55171-plywrite?s_tid=mwa_osa_a
function plyWrite(xyz,face,filename)
fid = fopen(filename,'w+');
fileBase = strsplit(filename,'.');
if min(min(face)) > 0
    face = face - min(min(face));
end
%% Write the Header
fprintf(fid,'%s\n','ply');
fprintf(fid,'%s\n','format ascii 1.0');
fprintf(fid,'%s\n','comment VTK generated PLY File');
fprintf(fid,'%s\n','obj_info vtkPolyData points and polygons: vtk4.0');
fprintf(fid,'%s%d\n','element vertex ',size(xyz,1));
fprintf(fid,'%s\n','property float x');
fprintf(fid,'%s\n','property float y');
fprintf(fid,'%s\n','property float z');
fprintf(fid,'%s%d\n','element face ',size(face,1));
fprintf(fid,'%s\n','property list uchar int vertex_indices');
fprintf(fid,'%s\n','end_header');
% for ii = 1:size(xyz,1)
%     fprintf(fid,'%f %f %f\n',xyz(ii,1),xyz(ii,2),xyz(ii,3));
% end
fprintf(fid,'%.15G %.15G %.15G\n',[xyz(:,1)';xyz(:,2)';xyz(:,3)']);
% for ii = 1:size(face,1)
%     fprintf(fid,'%d %d %d %d\n',3,face(ii,1),face(ii,2),face(ii,3));
% end
fprintf(fid,'%d %d %d %d\n',[(ones(size(face,1),1)*3)';face(:,1)';face(:,2)';face(:,3)']);
fclose(fid);
