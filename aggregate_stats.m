% repeated measures change in power across regions

% run aggregate power calc first
run pow_aggregateSpectral.m

% get change in power from baseline

nf_ix = [2:8];
nf_r = 4:6;
nf_n = [1:3,7];

for gr_ix = 1:length(groups)
    % p-val of absolute theta during NF
%     pval_abs(gr_ix) = friedman(thetaPW.(groups{gr_ix})(:,2:8),1);
    
    % change in theta from baseline
    thetaChange_agg.(groups{gr_ix}) = thetaPW.(groups{gr_ix}) ./ thetaPW.(groups{gr_ix})(:,1) - 1;
    thetaChange_agg.(groups{gr_ix})(:,1) = [];
    pval_change(gr_ix) = friedman(thetaChange_agg.(groups{gr_ix})(:,1:7),1);

    % pre-avgNF-post friedman
    alphaAvg_agg.(groups{gr_ix}) = [alphaPW.(groups{gr_ix})(:,1) mean(alphaPW.(groups{gr_ix})(:,2:8),2) alphaPW.(groups{gr_ix})(:,9)];
    thetaAvg_agg.(groups{gr_ix}) = [thetaPW.(groups{gr_ix})(:,1) mean(thetaPW.(groups{gr_ix})(:,2:8),2) thetaPW.(groups{gr_ix})(:,9)];
%     pval_fri_avg(gr_ix) = friedman(thetaAvg_agg.(groups{gr_ix}),1);

    % alpha-theta ratio
    alphaThetaRatio.(groups{gr_ix}) = 100*alphaPW.(groups{gr_ix}) ./ thetaPW.(groups{gr_ix});
%     pval_ratio(gr_ix) = friedman(alphaThetaRatio.(groups{gr_ix})(:,2:8),1);

    % alpha-theta ratio on avg NF
    ratioAvg_agg.(groups{gr_ix}) = [alphaThetaRatio.(groups{gr_ix})(:,1) mean(alphaThetaRatio.(groups{gr_ix})(:,2:8),2) alphaThetaRatio.(groups{gr_ix})(:,9)];
%     pval_ratioavg(gr_ix) = friedman(ratioAvg_agg.(groups{gr_ix}),1);

    [ath_corr.(groups{gr_ix})(1), ath_corr.(groups{gr_ix})(2)] = ...
        corr(thetaAvg_agg.(groups{gr_ix})(:,2), alphaAvg_agg.(groups{gr_ix})(:,2), 'type', 'Spearman');
    
        % corr of theta and alpha by subject
    for subidx = 1:length(thetaPW.(groups{gr_ix}))
        [ath_bySub.(groups{gr_ix})(subidx,1), ath_bySub.(groups{gr_ix})(subidx,2)] = ...
            corr(thetaPW.(groups{gr_ix})(subidx,:)', alphaPW.(groups{gr_ix})(subidx,:)', 'type', 'Spearman');
    end
end