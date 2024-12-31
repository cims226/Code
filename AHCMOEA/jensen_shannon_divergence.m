function jsd = jensen_shannon_divergence(P, Q)
%-------------------------------------------------------------------------%   
    
    P = P / sum(P);
    Q = Q / sum(Q);
    M = 0.5 * (P + Q); 
    KLD_PM = kldivergence(P, M);
    KLD_QM = kldivergence(Q, M); 
    jsd = 0.5 * (KLD_PM + KLD_QM);
end

%-------------------------------------------------------------------------%

function kld = kldivergence(P, Q)
    kld = sum(P .* log(P ./ Q), 'omitnan');
end



