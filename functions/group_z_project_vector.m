function [vector_zproj,std_zproj,SE_zproj] = group_z_project_vector(vector,group_z_value);

sz = size(vector);
sz = squeeze(sz);
vector = squeeze(vector);
slice_number = floor(sz/group_z_value);  % divides the total number of samples by the downsampling factor (in our case 500Hz to downsample to 1Hz)
slice_number = squeeze(slice_number);
%groupZ_value = sz(1)/slice_number;
vector2 = vector(1:slice_number*group_z_value);  %take a 'slice' of our data (really it's just taking the whole lot of data)
vector2 = reshape(vector2,[group_z_value,slice_number(1)]); %reshape vector2 into downsampling factor by slice_number
vector_zproj = mean(vector2,1); % takes the mean for every slice (i.e. every 500 samples of our original data)%
%SE_zproj = (std(vector2,[],3))./slice_number;  % Not sure- standard error for the same point over every slice?
vector_zproj = squeeze(vector_zproj); %Squeeze removes singleton dimensions
%std_zproj = squeeze(std_zproj);
%SE_zproj = squeeze(SE_zproj);
