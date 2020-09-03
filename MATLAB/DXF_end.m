function [ ] = DXF_end( fid )
%DXF_END Summary of this function goes here
%   Detailed explanation goes here
% fprintf(fid,'0\r\nENDSEC\r\n0\r\nEOF\r\n');
fid_j=fopen('TIAL.DXF','r');
while feof(fid_j)==0
	fprintf(fid,[fgetl(fid_j) '\r\n']);
end
fclose(fid_j);
fclose(fid);
end

