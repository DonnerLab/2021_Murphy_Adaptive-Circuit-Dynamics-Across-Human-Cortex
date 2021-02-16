
function [LPRout,surprise,scaled_prior] = accPsi_fast(LLRin,L_nm1,psi_n,Stype,pIn,H)

% Initializing LPR and surprise vectors
LPRout = LLRin(:,1) + interp1(L_nm1,psi_n,zeros(size(LLRin,1),1),'spline');
surprise = nan(size(LPRout));
scaled_prior = zeros(size(LPRout));

% Run belief updating
for s = 2:size(LLRin,2)
    scaled_prior(:,end+1) = nan(size(scaled_prior,1),1);
    scaled_prior(LPRout(:,end)>=min(L_nm1) & LPRout(:,end)<=max(L_nm1),end) = interp1(L_nm1,psi_n,LPRout(LPRout(:,end)>=min(L_nm1) & LPRout(:,end)<=max(L_nm1),end),'spline');
    scaled_prior(LPRout(:,end)<min(L_nm1),end) = psi_n(1);   % setting any posterior beliefs below lowest fitted phi to current parameter estimate for lowest fitted phi
    scaled_prior(LPRout(:,end)>max(L_nm1),end) = psi_n(end);   % setting any posterior beliefs above highest fitted phi to current parameter estimate for highest fitted phi
    
    LPRout(:,end+1) = LLRin(:,s)+scaled_prior(:,end);
end

scaled_prior(:,end+1) = nan(size(scaled_prior,1),1);
scaled_prior(LPRout(:,end)>=min(L_nm1) & LPRout(:,end)<=max(L_nm1),end) = interp1(L_nm1,psi_n,LPRout(LPRout(:,end)>=min(L_nm1) & LPRout(:,end)<=max(L_nm1),end),'spline');
scaled_prior(LPRout(:,end)<min(L_nm1),end) = psi_n(1);   % setting any posterior beliefs below lowest fitted phi to current parameter estimate for lowest fitted phi
scaled_prior(LPRout(:,end)>max(L_nm1),end) = psi_n(end);   % setting any posterior beliefs above highest fitted phi to current parameter estimate for highest fitted phi

% Calculate sample-wise surprise
if strcmp(Stype,'pCP')  % analytically-derived change-point probability
    pR = l2p([zeros(size(LPRout,1),1) LPRout],'p');
    pL = l2p([zeros(size(LPRout,1),1) LPRout],'n');
    
    for s = 1:size(LLRin,2)
        surprise(:,s) = (H.* ((squeeze(pIn(:,s,1)).*pL(:,s)) + (squeeze(pIn(:,s,2)).*pR(:,s)))) ./ ...
            ((H.* ((squeeze(pIn(:,s,1)).*pL(:,s)) + (squeeze(pIn(:,s,2)).*pR(:,s)))) + ((1-H).*((squeeze(pIn(:,s,2)).*pL(:,s)) + (squeeze(pIn(:,s,1)).*pR(:,s)))));
    end
    
elseif strcmp(Stype,'DY')  % Dayan & Yu-type surprise
    for s = 2:size(LLRin,2)
        surprise(sign(LPRout(:,s-1))==1,s) = (1./(exp(LPRout(sign(LPRout(:,s-1))==1,s))+1))./(1./(exp(LPRout(sign(LPRout(:,s-1))==1,s-1))+1));  % if belief before this sample was 'positive' (p1>p2)
        surprise(sign(LPRout(:,s-1))==-1,s) = (exp(LPRout(sign(LPRout(:,s-1))==-1,s))./(exp(LPRout(sign(LPRout(:,s-1))==-1,s))+1))...
                                                ./(exp(LPRout(sign(LPRout(:,s-1))==-1,s-1))./(exp(LPRout(sign(LPRout(:,s-1))==-1,s-1))+1));  % if belief before this sample was 'negative' (p2>p1)
        surprise(LPRout(:,s-1)==0,s) = nan;
    end
    surprise = log(surprise);
    
elseif strcmp(Stype,'signedL') % Surprise as signed difference in belief (not scaled by prior uncertainty)
    for s = 2:size(LLRin,2)
        surprise(sign(LPRout(:,s-1))==1,s) = LPRout(sign(LPRout(:,s-1))==1,s-1)-LPRout(sign(LPRout(:,s-1))==1,s);
        surprise(sign(LPRout(:,s-1))==-1,s) = LPRout(sign(LPRout(:,s-1))==-1,s)-LPRout(sign(LPRout(:,s-1))==-1,s-1);
        surprise(LPRout(:,s-1)==0,s) = nan;
    end
    
elseif strcmp(Stype,'absL') % Surprise as unsigned difference in belief (conceptually similar to KL divergence)
    for s = 1:size(LLRin,2)
        surprise(:,s) = abs(scaled_prior(:,s+1)-scaled_prior(:,s));
    end
    
elseif strcmp(Stype,'scaled_prior')  % Dayan & Yu-type surprise
    for s = 2:size(LLRin,2)
        surprise(sign(scaled_prior(:,s))==1,s) = (1./(exp(LPRout(sign(scaled_prior(:,s))==1,s))+1))./(1./(exp(scaled_prior(sign(scaled_prior(:,s))==1,s))+1));  % if belief before this sample was 'positive' (p1>p2)
        surprise(sign(scaled_prior(:,s))==-1,s) = (exp(LPRout(sign(scaled_prior(:,s))==-1,s))./(exp(LPRout(sign(scaled_prior(:,s))==-1,s))+1))...
                                                ./(exp(scaled_prior(sign(scaled_prior(:,s))==-1,s))./(exp(scaled_prior(sign(scaled_prior(:,s))==-1,s))+1));  % if belief before this sample was 'negative' (p2>p1)
        surprise(scaled_prior(:,s)==0,s) = nan;
    end
    surprise = log(surprise);
    
elseif strcmp(Stype,'conditional')  % surprise conditional on both prior belief and stimulus probability
    for s = 1:size(LLRin,2)
        surprise(:,s) = -log((squeeze(pIn(:,s,2)).*(1./(exp(scaled_prior(:,s))+1))) + ...
            (squeeze(pIn(:,s,1)).*(1-(1./(exp(scaled_prior(:,s))+1)))));
    end
    
elseif strcmp(Stype,'DY_prior_weighted')  % Prior-weighted, choice-combined Dayan & Yu-type surprise
    for s = 1:size(LLRin,2)
        surprise(:,s) = log(l2p(LPRout(:,s),'p')./l2p(scaled_prior(:,s),'p')).*log(l2p(scaled_prior(:,s),'n')./l2p(scaled_prior(:,s),'p')) + ...
                        log(l2p(LPRout(:,s),'n')./l2p(scaled_prior(:,s),'n')).*log(l2p(scaled_prior(:,s),'p')./l2p(scaled_prior(:,s),'n'));
    end
    
end