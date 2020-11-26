%% scripts for retreiving and calculating monthly GMWL anomalies

%% Retrieve all file names of the downloaded zip files
in_dir = 'E:\WRRTest\GMW_SK\GWLObs';
out_dir = 'E:\WRRTest\GMW_SK\Output_Daily';
gtws_file = 'E:\WRRTest\GRACE_JPL_Mascon_RL06\Version10\GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01CRIv01.nc';

%% read the GRACE TWS data from NETCDF file to get the time stamps
%ncdisp(gtws_file)
day = ncread(gtws_file,'time');
day_bounds = ncread(gtws_file,'time_bounds');

total_month = 158;

%% get file names 
zipfnames = dir(in_dir);
% folders = {zipfnames([zipfnames.isdir]).folder};
fnames = {zipfnames(~[zipfnames.isdir]).name};
nf = length(fnames);

%% read in all data
% [dt,val,id] = readvars(in_daily_file);
% unique_id = unique(id);
% nf = length(uniqe_id);

%% process data well bt well
for i =1:nf
    csvfn = [in_dir '\' fnames{i}];
    %% read in the data
    [dt,val] = readvars(csvfn);
    %% convert the time stamps to days
    dt_day = convertDateTimeToDays(dt);
    %% initialize output data 
    monthly_data = ones(total_month, 3)*(-32760);
    valid_month = 0;
    for j=1:total_month         
        %% determine the day ranges for the monthly dataset i
        day_bounds_j = day_bounds(:,j);
        day_start = day_bounds_j(1);
        day_end = day_bounds_j(2);
        dayj = day(j);

        idx_days = find(dt_day>=day_start & dt_day<day_end);
        if isempty(idx_days)
            monthly_data(j,1) = day(j);
            % monthly_data(j,2) = -32760;
        else %% retreive the data and average them
            mdata = val(idx_days);
            monthly_data(j,1) = day(j);
            monthly_data(j,2) = mean(mdata);
            valid_month = valid_month + 1;
        end
    end
    %% deduct the baseline
    if valid_month == total_month
        baseline = mean(monthly_data(:,2));
        monthly_data(:,3) = monthly_data(:,2) - baseline;
    end
    %% save the results to file
    outfn = [out_dir '\total_' num2str(valid_month) '_month_' fnames{i}];
    writematrix(monthly_data,outfn);
end