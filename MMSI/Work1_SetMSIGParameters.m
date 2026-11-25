% Create native triangular MSIG.ForPar; metres, density kg/m^3, G=1.
clear; clc;
SourceSpacing=1; ObservationSpacing=1; Height=1;
R2=6371000; R1=R2-30000; Density=2670;
AbsTol=1e-16; RelTol=1e-16;

assert(SourceSpacing>0 && SourceSpacing<=30 && ObservationSpacing>0);
assert(abs(round(180/ObservationSpacing)*ObservationSpacing-180)<1e-10);
assert(R1>0 && R2>R1 && R2+Height>0 && Height~=0 && Height~=R1-R2);
nv=round(8*pi/(sqrt(3)*(SourceSpacing*pi/180)^2));
i=(0:nv-1)'; z=1-2*(i+0.5)/nv; angle=i*pi*(3-sqrt(5));
vertices=[sqrt(1-z.^2).*cos(angle),sqrt(1-z.^2).*sin(angle),z];
faces=convhull(vertices);
flip=sum(vertices(faces(:,1),:).*cross(vertices(faces(:,2),:),vertices(faces(:,3),:),2),2)<0;
faces(flip,:)=faces(flip,[1,3,2]);
rho=repmat(Density,nv,1);
[lon,lat]=meshgrid(-180:ObservationSpacing:180-ObservationSpacing,...
    -90+ObservationSpacing:ObservationSpacing:90-ObservationSpacing);
coordinates=[[lon(:);0;0],[lat(:);-90;90],repmat(R2+Height,numel(lon)+2,1)];

fid=fopen('MSIG.ForPar','wb','ieee-le'); assert(fid>=0,'Cannot open parameter file');
cleanup=onCleanup(@()fclose(fid));
fwrite(fid,uint8([double('TSPMSI1'),0]),'uint8');
fwrite(fid,[AbsTol,RelTol],'double'); fwrite(fid,[nv,size(faces,1),1,0],'int32');
fwrite(fid,vertices.','double');
fwrite(fid,repmat(R2,nv,1),'double'); fwrite(fid,repmat(R1,nv,1),'double');
fwrite(fid,rho,'double'); fwrite(fid,faces.'-1,'int32');
fwrite(fid,size(coordinates,1),'int32'); fwrite(fid,coordinates,'double');
assert(isempty(ferror(fid)),'Failed to write parameter file'); clear cleanup;
fprintf('MSIG.ForPar: %d vertices, %d triangles, %d observations; AbsTol=%.1e RelTol=%.1e\n',...
    nv,size(faces,1),size(coordinates,1),AbsTol,RelTol);
