addpath('C:\Users\2146112s\Documents\PhDina\DataAnalysis-preproc\common_fct\')
addpath('C:\Users\2146112s\Documents\PhDina\DataAnalysis-preproc\neuralPaths\')
%% correlate with what?
corr_hr = input('Press \n1 to get correlation with FOOOF-extracted alpha power change \n2 for HR - Pearsons coeff \n3 for HR corr by strategy \n4 for Hurst-HR \n');

switch corr_hr 
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
        load([folder_results 'BandPower_IAF.mat']);
        subcode_pow = erase({PW_all.code},'sub');
end

% load HR data
folder_ax = 'C:\Users\2146112s\OneDrive - University of Glasgow\PhDina\healthy_study\assessments\HR_MEP\';

load([folder_ax 'HR\' modType '\allHR_NF.mat'])
load([folder_ax 'HR\' modType '\allHR_post.mat'])
load([folder_ax 'HR\' modType '\allHR_pre.mat'])
load([folder_ax 'HR\' modType '\allM_NF.mat'])
load([folder_ax 'HR\' modType '\allM_post.mat'])
load([folder_ax 'HR\' modType '\allM_pre.mat'])

subcode_hr = fieldnames(allHR_pre);
[~, fixd_ord] = ismember(subcode_hr,subcode_pow);

allHR = allHR_NF; allHR.pre = allHR_pre;  allHR.post = allHR_post;
allM = allM_NF; allM.pre = allM_pre;  allM.post = allM_post;

% get normalised H-reflex
cond = {'pre','trial1','trial2','trial3','post'};
hr_norm = getHRnorm(allHR, allM, cond, subcode_hr);

hr_norm_chg = 100*(hr_norm ./ hr_norm(:,1)-1);
hr_norm_chg(:,1) = [];
%% correlate
switch corr_hr
    case 4
        load(['C:\Users\2146112s\OneDrive - University of Glasgow\PhDina\healthy_study\analysis\HurstExp\H_all_' modType '.mat'])
        
        % get Hurst change
        hpre = H_exp.H(fixd_ord,1);
        h_nfhr = H_exp.H(fixd_ord,4:6);
        h_nfchg = 100*(h_nfhr./hpre -1);
       
        for tr_ix = 1:3
            % looks like the last one is an outlier
            mdl = linregmld_plot(hr_norm_chg(:,tr_ix),h_nfhr(:,tr_ix),['HR vs avg Hurst NF' int2str(tr_ix+3)],1);
        end
        % bpow
%         bp_pre = bandpow.alpha.PreEO(:,28);
%         bp_nf = [bandpow.alpha.NF4(:,28) bandpow.alpha.NF5(:,28) bandpow.alpha.NF6(:,28);];
%         bp_nf = 100*(bp_nf./bp_pre-1);
%         bp_nf = bp_nf(fixd_ord,:);
%         
%         for tr_ix = 1:3
%             success_ix = find(bp_nf(:,tr_ix) >10);
%             % looks like the last one is an outlier
%             mdl = linregmld_plot(hr_norm_chg(success_ix,tr_ix),h_nfhr(success_ix,tr_ix),['HR vs avg Hurst successful NF' int2str(tr_ix+3)],1);
%         end
end

