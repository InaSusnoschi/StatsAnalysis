addpath('C:\Users\2146112s\Documents\PhDina\DataAnalysis-preproc\common_fct\')

%% correlate with what?
corr_afooof = input('Press \n1 to get correlation with FOOOF-extracted alpha power change \n2 for MEP - Pearsons coeff \n3 for MEP corr by strategy \n4 for Hurst-MEP \n');
switch corr_afooof 
    case 1
        % correlation between FOOOF alpha and MEP change
        modType = 'up';
        addpath('C:\Users\2146112s\Documents\PhDina\DataAnalysis-preproc\FOOOFanalysis\')
        run step2_loadFOOOFgroup.m
        % get indices of included part
        subcode_pow = erase({fooof_files.name},'fooof_results_');
        subcode_pow = erase(subcode_pow,'.json');
        sub_success = subcode_pow(success_ix);
        
        % get average FOOOF alpha change
        achg_avg_nf     = 100*(a_nf_avg'./a_base(:,1) - 1);
        achg_avg_nfall  = 100*(a_nfall_avg'./a_base(:,1) - 1);
        achg_post       = 100*(a_post(:,1)./a_base(:,1) - 1);
    case 2
        % run healthy stats corr
        run healthy_stats_correlations.m
        subcode_pow = erase({PW_all.code},'sub');
    case 3
        run healthy_stats_correlations.m
        subcode_pow = erase({PW_all.code},'sub');
        mdt = split(modType,'\'); modType = mdt{1};
        run fct_getPartStrategy.m
    case 4
        expGroup = input('Insert 1 for alpha-UP group, 2 for alpha-DOWN group, 3 for both: ');

        switch expGroup
            case 1
                modType = 'up';
            case 2
                modType = 'down';
            case 3
                modType = 'updown';
        end
        folder_results  = (['C:\Users\2146112s\OneDrive - University of Glasgow\PhDina\healthy_study\analysis\spectral\' modType '\session3\' ]);
        load([folder_results 'Power_allSub.mat']);
        subcode_pow = erase({PW_all.code},'sub');
end


% load MEP data
folder_ax = 'C:\Users\2146112s\OneDrive - University of Glasgow\PhDina\healthy_study\assessments\HR_MEP\';

load([folder_ax 'MEP\' modType '\allMEP_pre.mat'])
load([folder_ax 'MEP\' modType '\allMEP_post1.mat'])
load([folder_ax 'MEP\' modType '\allMEP_post2.mat'])

% get only participants I have MEP for
subcode_mep = fieldnames(allMEP_pre);
[~, fixd_ord] = ismember(subcode_mep,subcode_pow);
if corr_afooof == 3
    % select strategy here
    subin_mep = find(ismember(subcode_mep,sub_decr));
    % intersect indices (MEP-pow, MEP-strategy)
    success_ix = intersect(fixd_ord,subin_mep');
    corr_afooof = 2; % same as case 2 for further correlations
end

%% get MEP change
meanz = zeros(3,length(subcode_pow)); stdmepz = meanz; meanz_scaled = meanz;
    p=zeros(1,length(subcode_pow));
    for subidx = 1:length(subcode_pow)

        if isfield(allMEP_post2, subcode_pow{subidx})
            maxlen = max([length(allMEP_pre.(subcode_pow{subidx})),length(allMEP_post1.(subcode_pow{subidx})),length(allMEP_post2.(subcode_pow{subidx}))]);
            allmeps = nan(maxlen,3);
            allmeps(1:length(allMEP_pre.(subcode_pow{subidx})),1) = allMEP_pre.(subcode_pow{subidx})';
            allmeps(1:length(allMEP_post1.(subcode_pow{subidx})),2) = allMEP_post1.(subcode_pow{subidx})';
            allmeps(1:length(allMEP_post2.(subcode_pow{subidx})),3) = allMEP_post2.(subcode_pow{subidx})';
        elseif isfield(allMEP_post1, subcode_pow{subidx})
            maxlen = max([length(allMEP_pre.(subcode_pow{subidx})),length(allMEP_post1.(subcode_pow{subidx}))]);
            allmeps = nan(maxlen,3);
            allmeps(1:length(allMEP_pre.(subcode_pow{subidx})),1) = allMEP_pre.(subcode_pow{subidx})';
            allmeps(1:length(allMEP_post1.(subcode_pow{subidx})),2) = allMEP_post1.(subcode_pow{subidx})';
        else 
            allmeps = nan(1,3);
        end
        allmep_bysub(subidx).pre = allmeps(:,1);
        allmep_bysub(subidx).post1 = allmeps(:,2);
        allmep_bysub(subidx).post2 = allmeps(:,3);
        
        meanz(:,subidx) = mean(allmeps,1,'omitnan');
        stdmepz(:,subidx) = std(allmeps,[],1,'omitnan');
        
        meanz_scaled(:,subidx) = mean(fct_robust_scale(allmeps),1,'omitnan');
    end
    
post_v_pre(1,:) = 100*(meanz(2,:)./meanz(1,:) - 1);
post_v_pre(2,:) = 100*(meanz(3,:)./meanz(1,:) - 1);

%% correlate
switch corr_afooof
    case 1
        % Correlation with alpha power change (from FOOOF)
        mepchg = post_v_pre(1,:);
        
        if length(subcode_mep) == length(subcode_pow)
            obs2 = achg_avg_nf;
            mdl = linregmld_plot(mepchg',obs2,'MEP vs FOOOF-alpha (all)',1);
            
            % successful
            obs2 = achg_avg_nf(success_ix);
            mdl = linregmld_plot(mepchg(success_ix)',obs2,'MEP vs FOOOF-alpha (successful)',1);
            
            % NF7
            obs2 = cell2mat(a_change(:,8));
            mdl = linregmld_plot(mepchg',obs2,'MEP vs FOOOF-alpha (NF7)',1);
            
            % successful NF7
            mdl = linregmld_plot(mepchg(success_ix)',obs2(success_ix),'MEP vs FOOOF-alpha (NF7)',1);
            
            % post - can't use it for those with no peak post-NF
        end

    case 2
        % Correlation with Pearsons or Spearmans coefficients
        if statCalc == 3 || statCalc == 4
                % pearsons corr coeff for MEP participants
                aa = {pears.rho};
                aa( cellfun(@isempty, aa) ) = {0};
                pearr = cell2mat(aa);
                fixd_ord(cellfun(@isempty,{pears.rho})) = []; % for FOOOF case
                
                % get linear regression with everyone
                mdl = linregmld_plot(pearr(fixd_ord)',post_v_pre(1,fixd_ord)','MEP vs Pearsons r',1);
                grid on; xline(0,'r'); yline(0,'r');
                xlabel('Correlation coefficient'); ylabel('MEP change from Pre [%]')
                
                % get successful part only (unless strategy corr calc)
                if ~exist('success_ix','var')
%                     mdt = split(modType,'\');
                    nfix = fct_getBestSubsession(modType,10,'abs');
                    success_ix = find(~isnan(nfix));
                    success_ix = intersect(fixd_ord,success_ix'); %should stay the same for non-FOOOF
                end
                
                pearr_s = pearr(success_ix);
                mep_1_s = post_v_pre(1,success_ix);
                mep_2_s = post_v_pre(2,success_ix);
                
                mdl = linregmld_plot(pearr_s',mep_1_s(1,:)','MEP vs Pearsons r (successful)',1);
                grid on; xline(0,'r'); yline(0,'r');
                xlabel('Correlation coefficient'); ylabel('MEP change from Pre [%]')
                
                if exist('spear','var')
                    aa = {spear.rho};
                    aa( cellfun(@isempty, aa) ) = {0};
                    spearr = cell2mat(aa);
                    mdl = linregmld_plot(spearr',post_v_pre(1,:)','MEP vs Spearman`s r',1);
                    grid on; xline(0,'r'); yline(0,'r');
                    xlabel('Correlation coefficient'); ylabel('MEP change from Pre [%]')
                    spearr_s = spearr(success_ix);
                    mdl = linregmld_plot(spearr_s',mep_1_s(1,:)','MEP vs Spearman r (successful)',1);
                end
        elseif statCalc == 5
                % run healthy stats corr for pre v post:
                a_chg = 100*(abs_alpha_post./abs_alpha_pre - 1);
                a_incr = find(a_chg<-10);
                mep_incr = post_v_pre(1,a_incr);
                post_incr = a_chg(a_incr);
                mdl = linregmld_plot(mep_incr',post_incr,'MEP vs Pearsons r (successful)',1);
                
                %hebbian
                a_heb = find(a_chg(1:23)>10);
                a_heb2 = find(a_chg(24:end)<-10);
                a_heb = [a_heb ; 23+a_heb2];
                mdl = linregmld_plot(post_v_pre(1,a_heb)',a_chg(a_heb),'MEP vs IA change post-NF(Hebbian)',1);
                
                %non-hebbian
                a_nheb = find(a_chg(1:23)<-10);
                a_nheb2 = find(a_chg(24:end)>10);
                a_nheb = [a_nheb ; 23+a_nheb2];
                mdl = linregmld_plot(post_v_pre(1,a_nheb)',a_chg(a_nheb),'MEP vs IA change post-NF(non-Hebbian)',1);
        end
    case 4
        load(['C:\Users\2146112s\OneDrive - University of Glasgow\PhDina\healthy_study\analysis\HurstExp\H_all_' modType '.mat'])
               
        for tr_ix = 2:9
%             h_chg = 100*(H_exp.H(:,tr_ix)./H_exp.H(:,1)-1);
            h_chg = H_exp.H(:,tr_ix);
            mdl = linregmld_plot(h_chg(fixd_ord),post_v_pre(1,fixd_ord)',[modType ': MEP (1) vs Hurst change NF' int2str(tr_ix-1)],1);
            grid on; %xline(0,'r'); yline(0,'r');
            xlabel('Hurst exponent'); ylabel('MEP change from Pre [%]')
        end
       nfix = fct_getBestSubsession(modType,10,'abs');
       success_ix = find(~isnan(nfix));
       success_ix = intersect(fixd_ord,success_ix');
       % average change in H-exp (in successful participants) vs MEP
       avg_hurst = mean(H_exp.H(success_ix,2:8),2);
       mdl = linregmld_plot(avg_hurst,post_v_pre(1,success_ix)','MEP vs avg Hurst (success)',1);
       % best subsession H-exp vs MEP
       best_hurst=[];
       for i = 1:length(nfix)
           if ~isnan(nfix(i))
            best_hurst = [best_hurst H_exp.H(i,nfix(i)+1)];
           end
       end
       mdl = linregmld_plot(post_v_pre(1,success_ix)',best_hurst','MEP vs best Hurst (success)',1);
       grid on;

end
