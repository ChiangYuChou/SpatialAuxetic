function [ fid ] = dxf_header( filename )
%DXF_HEADER Summary of this function goes here
%   Detailed explanation goes here
fid=fopen(filename,'w');
% fprintf(fid,'0\r\nSECTION\r\n');
fid_j=fopen('HEAD.DXF','r');
while feof(fid_j)==0
	fprintf(fid,[fgetl(fid_j) '\r\n']);
end
fclose(fid_j);
end

