function [b] = DXF_faces( fid,xyz,f,b)
%DXF_lines is set for export dxf file for SAP2000
%   fid is the file identity number generate by dxf_header
%	xyz is the corners of the faces
%	lin is the index number of the corners of the faces
for a=1:length(f)
	fprintf(fid,'0\r\n3DFACE\r\n5\r\nB%d\r\n100\r\nAcDbEntity\r\n8\r\nSHELLS\r\n100\r\nAcDbFace\r\n',a+b);
	%create new face element in layer SHELLS
	fprintf(fid,'10\r\n%.8f\r\n20\r\n%.8f\r\n30\r\n%.8f\r\n',xyz(f{a}(1),:));
	%first coordinate triple - first corner
	fprintf(fid,'11\r\n%.8f\r\n21\r\n%.8f\r\n31\r\n%.8f\r\n',xyz(f{a}(2),:));
	%second coordinate triple - second corner
	fprintf(fid,'12\r\n%.8f\r\n22\r\n%.8f\r\n32\r\n%.8f\r\n',xyz(f{a}(3),:));
	%third
    if length(f{a})>3
        fprintf(fid,'13\r\n%.8f\r\n23\r\n%.8f\r\n33\r\n%.8f\r\n',xyz(f{a}(4),:));
        %forth
    else
        fprintf(fid,'13\r\n%.8f\r\n23\r\n%.8f\r\n33\r\n%.8f\r\n',xyz(f{a}(3),:));
    end
end
b=a+b;
end
