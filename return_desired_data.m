function [data_desired]=return_desired_data(data,dnum,desired_period)
% This scripts selects time series for a given temporal window
data = double(data);
yr_initial = min(desired_period);
yr_final = max(desired_period);
ind = find(dnum>=datenum(yr_initial,1,1) & dnum<=datenum(yr_final,12,31));
data_desired = data(ind,:);
