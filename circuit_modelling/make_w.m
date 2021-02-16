clear, close all

wplus_vals = [1.68];
savepath = 'C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Analysis\Modelling_attractor\wang2002\wang2002lisa\';

seedval = 1;
if ischar(seedval); seedval = str2num(seedval);end
seedval = uint32(seedval);
st = RandStream('mt19937ar','Seed',seedval);
RandStream.setGlobalStream(st);

N           = 2000;  %Ncells
f           = 0.15;  %fraction selective

Ne = .8 * N ;   %number excitatory
Ni = .2 * N ;   %number inhibitory

%indexing:
ind.Pool1 = 1:(Ne*f) ;
ind.Pool2 = (Ne*f+1):(Ne*f*2) ;
ind.Ns = (Ne*f*2+1):Ne ;
ind.Inh = (Ne+1):N ;

%%
for i = 1:length(wplus_vals)
    fprintf(1, '\nBuilding Weight Matrix... ') ;
    w = zeros(N) ;  %weights matrix - all to all
    wplus = wplus_vals(i);
    wminus = 1 - f * (wplus - 1)/(1 - f) ;
    
    for i=1:N
        for j=1:N
            if i==j % No Self Connections
                w(i,j) = 0 ;
            elseif any(i == ind.Pool1) && any(j == ind.Pool1) % Same group
                w(i,j) = wplus ;
            elseif any(i == ind.Pool2) && any(j == ind.Pool2) % Same group
                w(i,j) = wplus ;
            elseif any(i == ind.Pool1) && any(j == ind.Pool2)  %% Different Group
                w(i,j) = wminus ;
            elseif any(i == ind.Pool2) && any(j == ind.Pool1) %% Different Group
                w(i,j) = wminus ;
            elseif any(i == union(ind.Pool1,ind.Pool2)) && any(j == ind.Ns) %Nonselective
                w(i,j) = wminus ;
            else %all other connections
                w(i,j) = 1;
            end
        end
    end
    if 0 %this is how to create the sparse matrix
        indr = randsample(N^2,N^2/2,false); %randomly sample half the connections
        w(indr) = 0; %set these to 0
        for i = 1:N;
            w(:,i) = w(:,i)./mean(w(:,i)); %normalise so mean input is still 1
        end
    end
    save([savepath,'w',num2str(wplus),'.mat'],'w')
end
