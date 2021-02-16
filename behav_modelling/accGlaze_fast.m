% Creates sequence of updated beliefs given some input sequence of log
% likelihood ratios, according the Glaze et al. normative model. H is
% hazard rate of changes in generative distribution, startpoint is
% starting belief in bans (LLR units).

function [LPRout,surprise,scaled_prior] = accGlaze_fast(LLRin,H,startpoint,Stype,pIn)

% Initializing LPR and surprise vectors
if H>0
    LPRout = LLRin(:,1)+startpoint+log(((1-H)/H)+exp(-startpoint))-log(((1-H)/H)+exp(startpoint));
else
    LPRout = LLRin(:,1)+startpoint;
end
surprise = nan(size(LPRout));
scaled_prior = zeros(size(LPRout));

% Run belief updating
if H>0
    for s = 2:size(LLRin,2)
        scaled_prior(:,end+1) = LPRout(:,end)+log(((1-H)/H)+exp(-LPRout(:,end)))-log(((1-H)/H)+exp(LPRout(:,end)));
        LPRout(:,end+1) = LLRin(:,s)+scaled_prior(:,end);
    end
else
    for s = 2:size(LLRin,2)
        scaled_prior(:,end+1) = LPRout(:,end);
        LPRout(:,end+1) = LLRin(:,s)+LPRout(:,end);
    end
end

% Calculate sample-wise surprise
if strcmp(Stype,'pCP')  % analytically-derived change-point probability
    pR = l2p([zeros(size(LPRout,1),1) LPRout],'p');
    pL = l2p([zeros(size(LPRout,1),1) LPRout],'n');
    
    for s = 1:size(LLRin,2)
        surprise(:,s) = (H.* ((squeeze(pIn(:,s,1)).*pL(:,s)) + (squeeze(pIn(:,s,2)).*pR(:,s)))) ./ ...
            ((H.* ((squeeze(pIn(:,s,1)).*pL(:,s)) + (squeeze(pIn(:,s,2)).*pR(:,s)))) + ((1-H).*((squeeze(pIn(:,s,2)).*pL(:,s)) + (squeeze(pIn(:,s,1)).*pR(:,s)))));
    end
    
elseif strcmp(Stype,'pCPpsi')  % analytically-derived change-point probability, but using PSI rather than L
    pR = l2p(scaled_prior,'p');
    pL = l2p(scaled_prior,'n');
    
    if H<=0.5,
        Heff = H;  % making expected rate of change relative to my TRANSFORMED expectation
    else Heff = 1-H;
    end
    
    for s = 1:size(LLRin,2)
        surprise(:,s) = (Heff.* ((squeeze(pIn(:,s,1)).*pL(:,s)) + (squeeze(pIn(:,s,2)).*pR(:,s)))) ./ ...
            ((Heff.* ((squeeze(pIn(:,s,1)).*pL(:,s)) + (squeeze(pIn(:,s,2)).*pR(:,s)))) + ((1-Heff).*((squeeze(pIn(:,s,2)).*pL(:,s)) + (squeeze(pIn(:,s,1)).*pR(:,s)))));
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
    for s = 2:size(LLRin,2)
        surprise(:,s) = abs(LPRout(:,s)-LPRout(:,s-1));
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
        surprise(:,s) = log(l2p(LPRout(:,s),'p')./l2p(scaled_prior(:,s),'p')).*(log(l2p(scaled_prior(:,s),'n')./l2p(scaled_prior(:,s),'p'))) + ...
            log(l2p(LPRout(:,s),'n')./l2p(scaled_prior(:,s),'n')).*(log(l2p(scaled_prior(:,s),'p')./l2p(scaled_prior(:,s),'n')));
    end
    
end