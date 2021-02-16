function pout = l2p(Lin,dir)

if strcmp(dir,'n')
    pout = 1./(exp(Lin)+1);
elseif strcmp(dir,'p')
    pout = exp(Lin)./(exp(Lin)+1);
end