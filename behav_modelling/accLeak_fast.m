% Creates sequence of updated beliefs given some input sequence of log
% likelihood ratios, according the Glaze et al. normative model.

function [LPRout,surprise,prior] = accLeak_fast(LLRin,leak,startpoint,Stype,pIn,H)

% Initializing LPR and surprise vectors
LPRout = LLRin(:,1)+startpoint;
surprise = nan(size(LPRout));

% Run belief updating %
for s = 2:size(LLRin,2)
    LPRout(:,end+1) = LPRout(:,end).*(1-leak)+LLRin(:,s);
end
prior = [zeros(size(LPRout,1),1) LPRout(:,1:end-1)];

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
    for s = 2:size(LLRin,2)
        surprise(:,s) = abs(LPRout(:,s)-LPRout(:,s-1));
    end
    
elseif strcmp(Stype,'DY_prior_weighted')  % Prior-weighted, choice-combined Dayan & Yu-type surprise
    for s = 2:size(LLRin,2)
        surprise(:,s) = log(l2p(LPRout(:,s),'p')./l2p(LPRout(:,s-1),'p')).*(log(l2p(LPRout(:,s-1),'n')./l2p(LPRout(:,s-1),'p'))) + ...
            log(l2p(LPRout(:,s),'n')./l2p(LPRout(:,s-1),'n')).*(log(l2p(LPRout(:,s-1),'p')./l2p(LPRout(:,s-1),'n')));
    end
end
