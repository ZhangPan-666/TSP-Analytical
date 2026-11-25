function result=Work2_ReadGravityForwardResults(file)
% Raw G=1 gravity derivatives in the fixed observer north/east/up frame.
if nargin==0, file='ResultG.dat'; end
names={'V','Vx','Vy','Vz','Vxx','Vxy','Vyy','Vzx','Vzy','Vzz',...
    'Vxxx','Vxxy','Vxxz','Vxyz','Vyyx','Vyyy','Vyyz','Vzzx','Vzzy','Vzzz'};
fid=fopen(file,'rb','ieee-le'); assert(fid>=0,'Cannot open result file');
cleanup=onCleanup(@()fclose(fid));
n=fread(fid,1,'int32=>double'); assert(isscalar(n) && n>=0,'Invalid observation count');
fseek(fid,0,'eof'); assert(ftell(fid)==4+8*n*(3+numel(names)),'Invalid result file size');
fseek(fid,4,'bof'); data=reshape(fread(fid,n*(3+numel(names)),'double=>double'),n,3+numel(names));
result=struct('Longitude',data(:,1),'Latitude',data(:,2),'Radius',data(:,3));
for k=1:numel(names), result.(names{k})=data(:,k+3); end
end
