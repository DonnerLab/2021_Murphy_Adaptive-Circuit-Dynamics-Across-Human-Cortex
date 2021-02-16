function value = default_arguments(argument_list, name, default)
value = default;
for i = 1:length(argument_list)
    if strcmp(argument_list{i}, name)
        value = argument_list{i+1};
    end
end

