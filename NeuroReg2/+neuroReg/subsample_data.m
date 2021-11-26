function data_sub=subsample_data(data, Sampling_factor)
data_sub = data;
data_sub.value= imresize3(data.value, Sampling_factor);
data_sub.x = linspace(data.x(1,1)/Sampling_factor,data.x(1,end), size(data_sub.value,1)) ;
data_sub.y = linspace(data.y(1,1)/Sampling_factor,data.y(1,end), size(data_sub.value,2));

if size(data.value,3)>1
    data_sub.z = linspace(data.z(1,1)/Sampling_factor,data.z(1,end), size(data_sub.value,3));
end

end
