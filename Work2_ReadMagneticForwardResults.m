clc;
clear;
close all;

fid=fopen('Result_TFM.dat','rb');

NumCoordinates=fread(fid,1,'int32');

Longitude=fread(fid,NumCoordinates,'double');
Latitude=fread(fid,NumCoordinates,'double');
Radius=fread(fid,NumCoordinates,'double');

V=fread(fid,NumCoordinates,'double');

Vx=fread(fid,NumCoordinates,'double');
Vy=fread(fid,NumCoordinates,'double');
Vz=fread(fid,NumCoordinates,'double');

Vxx=fread(fid,NumCoordinates,'double');
Vxy=fread(fid,NumCoordinates,'double');
Vyy=fread(fid,NumCoordinates,'double');
Vzx=fread(fid,NumCoordinates,'double');
Vzy=fread(fid,NumCoordinates,'double');
Vzz=fread(fid,NumCoordinates,'double');

fclose(fid);

DataSave=[Longitude,Latitude,Radius,V,Vx,Vy,Vz,...
    Vxx,Vxy,Vyy,Vzx,Vzy,Vzz];
writematrix(DataSave,'MagneticForwardResults.dat');
