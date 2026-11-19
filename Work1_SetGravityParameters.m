clc; clear; close all;
% Each row is one spherical triangle. Vertices use latitude/longitude degrees.
% Boundaries are shortest GREAT-CIRCLE arcs, not latitude small-circle arcs.
Fai1=[90-0.50*180/pi]; Fai2=[90-0.63*180/pi]; Fai3=[90-1.05*180/pi];
Lamda1=[0]; Lamda2=[0.85*180/pi]; Lamda3=[0.39*180/pi];
R1=[1]; R2=[2]; rho=[1];
Longitude=[-25;0]; Latitude=[61;90]; Radius=[3.4;3.0];
model={Fai1,Fai2,Fai3,Lamda1,Lamda2,Lamda3,R1,R2,rho};
M=numel(R1); N=numel(Radius);
assert(all(cellfun(@numel,model)==M),'Model arrays must have equal lengths');
assert(numel(Longitude)==N&&numel(Latitude)==N,'Coordinate arrays must have equal lengths');
fid=fopen('TFG.ForPar','wb','ieee-le'); assert(fid>=0,'Cannot open parameter file');
cleanup=onCleanup(@()fclose(fid));
fwrite(fid,M,'int32'); for k=1:numel(model), fwrite(fid,model{k}(:),'double'); end
fwrite(fid,N,'int32'); fwrite(fid,Longitude(:),'double'); fwrite(fid,Latitude(:),'double'); fwrite(fid,Radius(:),'double');
clear cleanup;
