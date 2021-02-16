function spikes = spikeloop(dt, tstop, Ne, Ni, w, Stim)

N= Ne + Ni ;
VL = -70 ;
Vth = -50 ;
VE = 0 ;
VI = -70 ;
Vreset = -55 ;
Cm = 1000 * [0.5 * ones(Ne,1); 0.2 * ones(Ni,1)] ;
gL = [25 * ones(Ne,1); 20 *  ones(Ni,1)] ;
tauref = [2 * ones(Ne,1); 1 * ones(Ni,1)] ;
taum = Cm ./ gL ;
U = 0.15; %LH
tauF = 2000; %LH
u = ones(Ne,1).*U;

gextAMPA = [2.1 * ones(Ne,1); 1.62 * ones(Ni,1)] ;
grecAMPA = [104/N * ones(Ne,1); 81/N * ones(Ni,1)] ;
gNMDA = [327/N * ones(Ne,1); 258/N * ones(Ni,1)] ;
gGABA = [1.3 * ones(Ne,1); 1 * ones(Ni,1)] ;
tauAMPA = 2 ;
tauNMDArise = 2 ;
tauNMDAdecay = 100 ;
tauGABA = 5 ;
alpha = 0.5  ;
Mg = 1 ;

Iw = w(:,(Ne+1):N);
Ew = w(:,1:Ne) ;
spikes = zeros(N,tstop) ;
Vt = Vreset * ones(N,1) ;
Isyn = zeros(N,1) ;
IextAMPA = zeros(N,1) ;
IrecAMPA = zeros(N,1) ;
IrecNMDA = zeros(N,1) ;
IrecGABA = zeros(N,1) ;
sAMPA = zeros(Ne,1) ;
sNMDA = zeros(Ne,1) ;
x = zeros(Ne,1) ;
sGABA = zeros(Ni,1) ;
sAMPAext = zeros(N,1) ;

for t = 2:tstop
    spikes(:,t) = Vt >= Vth ;
    Vt(Vt >= Vth) = Vreset ;
        
    sGABA  = sGABA .* exp(-dt/tauGABA) + spikes((Ne+1):N,t) ;     
    x = x .* exp(-dt/tauNMDArise) + spikes(1:Ne,t) ;
    
    sNMDA = sNMDA + (-sNMDA/tauNMDAdecay + alpha * x .* (1 - sNMDA)) * dt ;   
    sAMPA = sAMPA * exp(-dt/tauAMPA) + spikes(1:Ne,t) ;
    sAMPAext = sAMPAext * exp(-dt/tauAMPA) + Stim(:,t);
    u = u + dt*(U-u)./tauF + U*(1-u).*spikes(1:Ne,t); %LH
    ucache(:,t) = u;
    
    IextAMPA = gextAMPA .* (Vt - VE) .* sAMPAext ;
    IrecGABA = gGABA .* (Vt - VI) .* (Iw * sGABA) ;   
    IrecNMDA = (gNMDA .* (Vt - VE) ./ (1 + Mg / 3.57 * exp(1).^(-0.062*Vt))) .* (Ew * sNMDA) ;  
    IrecAMPA = grecAMPA .* (Vt - VE) .* (Ew * sAMPA) ;
    
    Isyn = (IextAMPA + IrecAMPA + IrecNMDA + IrecGABA)  ; 

    dV = (-gL .* (Vt - VL) - Isyn) ./ Cm   ;    
    Vt = Vt + dV * dt ;
        
    Vt(sum(spikes(:,max(t-floor(tauref/dt),1):t),2)>0) = Vreset ;
end


end