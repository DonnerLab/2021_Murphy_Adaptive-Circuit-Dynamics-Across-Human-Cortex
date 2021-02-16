function xout = dist_scatters(y,min_y)

% min_y: minimum allowed distance between dots on y dimension - if smaller than this, then offset along x

max_n_x = 6;   % maximum number of offset dots aong x dimension - if exceeds, then layer

y = y-min(y); y = y./max(y);  % normalize y to bounded b/w 0,1

[y,y_i] = sort(y);

x=[];
i=0;
while i<length(y)  % loop through each data point
	i = i+1;
    noff = length(find(y>=y(i) & y<=(y(i)+min_y))); % get number of dots to offset for this round
    
    ngrps = ceil(noff/max_n_x);  % number of sub-groups for this round
    if ngrps == 1  % if there's only 1 sub-group
        if noff == 1
            cx = 0.5;
        else cx = linspace(-0.5,0.5,noff).*(noff/max_n_x);
            [~,cxi] = sort(abs(cx));  % ensuring smallest y value will be most centered on x-axis
            cx = cx(cxi)+0.5;
        end
        x = [x cx];
    else  % if there's more than 1 sub-group
        for g = 1:ngrps
            if g == ngrps
                cnoff = mod(noff,max_n_x);
                if cnoff == 1
                    cx = 0.5;
                else cx = linspace(-0.5,0.5,cnoff).*(cnoff/max_n_x);
                    [~,cxi] = sort(abs(cx));
                    cx = cx(cxi)+0.5;
                end
            else cx = linspace(-0.5,0.5,max_n_x);
                [~,cxi] = sort(abs(cx));
                cx = cx(cxi)+0.5;
            end
            x = [x cx];
        end
    end
    i = i+noff-1;
end

xout=[];
for k = 1:length(x)
    xout(k,1) = x(y_i==k); % sorting x positions to align with input y
end
    
