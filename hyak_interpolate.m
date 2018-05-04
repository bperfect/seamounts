%% Convert variables in terrain following coordinates to gridded variables
%% Uses parallel computing toolbox for speed
%% 15 Jan 18
%% has issues closing the file. Need to close matlab to force the file to close

%% Obtain flow field in physical coordinates
clear
readAll = 1;
interpolating = 1;
directory = '/gscratch/stf/bperfect/ini_files/';
% directory = '/home/bperfect/seamounts/ROMS/seamount/';
% directory = '/media/bperfect/91d7719f-88e9-4c43-8698-be629b97a100/f0n1e3/';
%directory = 'C:\Users\bperfect\Documents\UW\Seamounts\ROMS\seamount_data\f1e5n1e3\';
casename = 'Case simname';
refinement = 3;

runname = 'simname';
hisname = [directory 'ocean_his_' num2str(refinement) 'x.nc'];
avgname = [directory 'ocean_avg_' num2str(refinement) 'x.nc'];
hisname = [directory 'ocean_his_' runname '.nc'];
avgname = [directory 'ocean_avg_' runname '.nc'];
%hisname = [directory 'ocean_his_f1e4n1e5.nc'];
%avgname = [directory 'ocean_avg_f1e4n1e5.nc'];


fh=fopen('simname','at');
fprintf(fh, 'Beginning interpolation to physical coordinates\n');

%% Parameters for sampling data
start=1; %time index to start sampling
nSamples = Inf;

if readAll == 1
    xx=1; yy=1; xstart=1; xsamples=Inf; ystart=1; ysamples=Inf;
else
    xx=1; yy=1; xstart=1; xsamples=Inf; ystart=1; ysamples=Inf;
end

dt = 1;
dz=50;
z_query = 0:-dz:-5000;
   
%% vertical grid variables
s_w         = ncread(hisname,'s_w');
s_rho       = ncread(hisname,'s_rho');
zeta        = 0;
hc          = ncread(hisname,'hc');
Cs_w        = ncread(hisname,'Cs_w');
Cs_r        = ncread(hisname,'Cs_r');

%% Horizontal grid variables
y_rho       = ncread(hisname,'y_rho',[xstart ystart],[xsamples ysamples],[xx yy]);
x_rho       = ncread(hisname,'x_rho',[xstart ystart],[xsamples ysamples],[xx yy]);
h           = ncread(hisname,'h',[xstart ystart],[xsamples ysamples],[xx yy]);
[Hz,z_w,z_r]= get_depth_Hz_ROMS(h,hc,s_rho,s_w,Cs_r,Cs_w,zeta(:,:,1),2);
z_r     = convert2rho(shiftdim(z_r,1),'rho');
z_w     = convert2rho(shiftdim(z_w,1),'rho');
W_r = interp1(squeeze(z_r(1,1,:)),eye(length(squeeze(z_r(1,1,:)))),z_query','spline',0);
W_w = interp1(squeeze(z_w(1,1,:)),eye(length(squeeze(z_w(1,1,:)))),z_query','spline',0);


%% Time variables
t=ncread(hisname,'ocean_time',start,nSamples,dt);
% t=t/24/3600; %convert from seconds to days

%% Step 1: obtain time means and delete full variables
fprintf(fh,'Loading variables for processing\n');

vars_to_interpolate = {'u','v','w','rho','AKt','AKv','tke','gls','rvorticity'};
varType = {'u','v','rho','rho','rho','rho','rho','rho','psi'};
filename = {hisname,hisname,hisname,hisname,hisname,hisname,hisname,hisname,avgname};
outputname = [directory '/interpolated_vars_' runname '.nc'];
% outputname = [directory 'interpolated_vars_tf.nc'];
%outputname = [directory 'interpolated_vars_.nc'];

terrain_type = {'r','r','w','r','w','w','w','w','r'};

%% Get the dimensions set up for writing netcdf files
x_rho = convert2rho(x_rho,'rho');
y_rho = convert2rho(y_rho,'rho');
h = convert2rho(h,'rho');
zeta = convert2rho(ncread(hisname,'zeta'),'rho');

dim1 = length(x_rho(:,1));
dim2 = length(x_rho(1,:));
dim3 = length(z_query);
dim4 = length(t);

ncid = netcdf.create(outputname,'NETCDF4')

% Define the dimensions of the variable.
dimid1 = netcdf.defDim(ncid,'x_rho',dim1);
dimid2 = netcdf.defDim(ncid,'y_rho',dim2);
dimid3 = netcdf.defDim(ncid,'z',dim3);
dimid4 = netcdf.defDim(ncid,'t',netcdf.getConstant('NC_UNLIMITED'));
dimid5 = netcdf.defDim(ncid,'timesteps',dim4);
varids = [];
for iteration = 1:length(vars_to_interpolate)
    varids = [varids netcdf.defVar(ncid,vars_to_interpolate{iteration},'NC_FLOAT',[dimid1,dimid2,dimid3,dimid4])];
    netcdf.defVarDeflate(ncid,varids(end),false,true,7);
end

varid1 = netcdf.defVar(ncid,'x_rho','NC_FLOAT',[dimid1,dimid2]);
varid2 = netcdf.defVar(ncid,'y_rho','NC_FLOAT',[dimid1,dimid2]);
varid3 = netcdf.defVar(ncid,'h','NC_FLOAT',[dimid1,dimid2]);
varid4 = netcdf.defVar(ncid,'t','NC_FLOAT',[dimid5]);
varid5 = netcdf.defVar(ncid,'z_query','NC_FLOAT',[dimid3]);
varid6 = netcdf.defVar(ncid,'zeta','NC_FLOAT',[dimid1,dimid2,dimid5]);

netcdf.endDef(ncid);

netcdf.putVar(ncid,varid1,x_rho);
netcdf.putVar(ncid,varid2,y_rho);
netcdf.putVar(ncid,varid3,h);
netcdf.putVar(ncid,varid4,t);
netcdf.putVar(ncid,varid5,z_query');
netcdf.putVar(ncid,varid6,zeta);
ncid

for iteration = 1:length(vars_to_interpolate)
%     fprintf(fh,['Interpolating ' vars_to_interpolate{iteration} ':\n']);
%     fprintf(fh,'[=========100========]\n');
%     fprintf(fh,'[');
  %  tic
    if terrain_type{iteration} == 'r'
        W = W_r;
        z = z_r;
    else
        W=W_w;
        z=z_w;
    end

    grid_var = zeros(dim1,dim2,dim3);
    fivePercent = floor(0.05*length(t));
    
    % rvorticity contains one fewer time index than everything else, so we
    % need to account for it
    if strcmp(vars_to_interpolate{iteration},'rvorticity')
        maxT = length(t)-1;
    else
        maxT = length(t);
    end
    
    for time = 1:maxT
        terrain_var = ncread(filename{iteration},vars_to_interpolate{iteration},[1 1 1 time],[Inf Inf Inf 1],[1 1 1 1]);
        terrain_var = convert2rho(terrain_var,varType{iteration});
        if mod(time,fivePercent) == 0
%             fprintf(fh,'=');
            disp(['5% checkpoint for ' vars_to_interpolate{iteration}])
        end
        %disp('time iteration');
        parfor (row = 1:dim1,28)
            for col = 1:dim2
                grid_var(row,col,:) = interp1(squeeze(z(row,col,:)),squeeze(terrain_var(row,col,:)),z_query','spline',0);
            end
        end
        netcdf.putVar(ncid,varids(iteration),[0 0 0 time],[dim1 dim2 dim3 1],grid_var);
        if mod(time,100)==0
            time
        end
    end
    % indicate that this variable is done and re-enter define mode
%     fprintf(fh,']... Complete!\n');
    
  %  toc
end
ncid
exit
netcdf.close(ncid);
fclose(fh);


