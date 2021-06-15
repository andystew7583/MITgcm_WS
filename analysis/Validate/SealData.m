%%%%%%%%%% Seal data

%%%%%%% 3 relevant datasets are ct70,ct43,ct27

datadir = 'data1/MITgcm_WS/analysis/Validate/UKSealData/ncARGO/';



for num = 574:866
    filename = fullfile(datadir,['ct43-' num2str(num) '-09_prof.nc']);
    % create a cell array of variables to load
    variables_to_load = {'LATITUDE','LONGITUDE','PRES_ADJUSTED','TEMP','PSAL_ADJUSTED'};
    if exist(filename)==2
% loop over the variables
        for j=1:numel(variables_to_load)
    % extract the jth variable (type = string)
            var = variables_to_load{j};

    % use dynamic field name to add this to the structure
            ct_43.(var) = ncread(filename,var);

    % convert from single to double
            if isa(ct_43.(var),'single')
                ct_43.(var) = double(ct_43.(var));
            end
        end
    end
end


for num = 1:5
    filename = fullfile(datadir,['ct27-W' num2str(num) '-07_prof.nc']);
    % create a cell array of variables to load
    variables_to_load = {'LATITUDE','LONGITUDE','PRES_ADJUSTED','TEMP','PSAL_ADJUSTED'};
    if exist(filename)==2
% loop over the variables
        for j=1:numel(variables_to_load)
    % extract the jth variable (type = string)
            var = variables_to_load{j};

    % use dynamic field name to add this to the structure
            ct_27.(var) = ncread(filename,var);

    % convert from single to double
            if isa(ct_27.(var),'single')
                ct_27.(var) = double(ct_27.(var));
            end
        end
    end
end


%%%%%%%%%% trying to interpolate seal data?

%%%%%% First have to load experiment data





run ../../newexp/defineGrid.m
run ../setExpname;
run ../loadexp;

dumpFreq = abs(diag_frequency(7));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);


%%%%%% Loading our Temp/Salinity data from respective experiment

tmin = 8*86400*365;
tmax = 9*86400*365;
Theta = readIters(exppath,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,1);
Salt = readIters(exppath,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,1);



N_Stations = 161;
ct_43LAT = ct_43.LATITUDE;
ct_43LON = ct_43.LONGITUDE;


temp_seal_c27 = zeros(length(zz),size(ct_43.PRES_ADJUSTED,1));
for n=1:size(ct_43.PRES_ADJUSTED,2)
  idx = find(~isnan(ct_43.PRES_ADJUSTED(:,n)));
  [b,i,j]=unique(ct_43.PRES_ADJUSTED(idx,n));
  temp_seal_c27(:,n) = interp1((b),ct_43.TEMP(i,n),squeeze(-zz),'linear');
end

Temp_new_lon = zeros(N_Stations,length(ymc),length(zz));
for k = 1:Ny
    for n=1:length(zz)
        seal_t = Theta(:,k,n);
        Temp_new_lon(:,k,n) = interp1(xmc,seal_t,ct_43.LONGITUDE,'linear');
    end
end

Temp_MOD = zeros(N_Stations,N_Stations,length(zz));
for k = 1:N_Stations
    for n=1:length(zz)
        seal_t = Temp_new_lon(k,:,n);
        Temp_MOD(k,:,n) = interp1(ymc,seal_t,ct_43.LATITUDE,'linear');
    end
end

hFacC_new_lon = zeros(N_Stations,length(ymc),length(zz));
for k = 1:Ny
    for n=1:length(zz)
        seal_hFacC = hFacC(:,k,n);
        hFacC_new_lon(:,k,n) = interp1(xmc,seal_hFacC,ct_43.LONGITUDE,'linear');
    end
end

hFacC_MOD = zeros(N_Stations,N_Stations,length(zz));
for k = 1:N_Stations
    for n=1:length(zz)
        seal_hFacC_lat = hFacC_new_lon(k,:,n);
        hFacC_MOD(k,:,n) = interp1(ymc,seal_hFacC_lat,ct_43.LATITUDE,'linear');
    end
end



[Seal_LO,MZ] = meshgrid(ct_43.LATITUDE,zz);

Temp_MOD = squeeze(Temp_MOD(end,:,:));

ylim([500 0])



