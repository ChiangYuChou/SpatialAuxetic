function [b] = DXF_lines( fid,xyz,lin,b)
%DXF_lines is set for export dxf file for SAP2000
%   fid is the file identity number generate by dxf_header
%	xyz is the points of the end of lines
%	lin is the index number of the two ends of each line
for a=1:length(lin(:,1))
	fprintf(fid,'0\r\nLINE\r\n5\r\nB%d\r\n100\r\nAcDbEntity\r\n8\r\nFRAMES\r\n100\r\nAcDbLine\r\n',a+b);
	%create new line element in layer FRAMES
	fprintf(fid,'10\r\n%.8f\r\n20\r\n%.8f\r\n30\r\n%.8f\r\n',xyz(lin(a,1),:));
	%first coordinate triple - starting point of line-element
	fprintf(fid,'11\r\n%.8f\r\n21\r\n%.8f\r\n31\r\n%.8f\r\n',xyz(lin(a,2),:));
	%second coordinate triple - ending point of line-element
end
b=a+b;
end

