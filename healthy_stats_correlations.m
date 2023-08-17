addpath('C:\Users\2146112s\Documents\PhDina\DataAnalysis-preproc\updated21')
addpath('C:\Users\2146112s\Documents\PhDina\DataAnalysis-preproc\common_fct')
% run Z_selectDataToProcess.m
expGroup = input('Insert 1 for alpha-UP group, 2 for alpha-DOWN group, 3 for both: ');

switch expGroup
    case 1
        modType = 'up';
    case 2
        modType = 'down';
    case 3
        modType = 'updown';
    otherwise
        modType = [];
end
% get successful participants
nfix = fct_getBestSubsession(modType,10,'abs');
success_ix = find(~isnan(nfix));

s_no = 3; % only calculate on session 3 data for now
chix = 28; % Cz

folder_results  = (['C:\Users\2146112s\OneDrive - University of Glasgow\PhDina\healthy_study\analysis\spectral\' modType '\session' int2str(s_no) filesep]);

load([folder_results 'BandPower_IAF.mat']);
% load([folder_results 'RelBandPower_IAF.mat']);
load([folder_results 'Power_allSub.mat']);

Subj_code = {PW_all.code};
% define trials etc.
all_trials =fieldnames(bandpow.alpha);
% trials_nf = {'NF1','NF2','NF3','NF4','NF5','NF6','NF7'};

trls = input(['Press \n1 to include baseline and all trials' ...
            '\n2 to include baseline and 1-3' ...
            '\n3 to include only NF1-3 \n4 Pre - Post']);
switch trls
    case 1
        trials_nf = {'PreEO','NF1','NF2','NF3','NF4','NF5','NF6','NF7'};
    case 2
        trials_nf = {'PreEO','NF1','NF2','NF3'};
    case 3
        trials_nf = {'NF1','NF2','NF3'};
    case 4
        trials_nf = {'PreEO','PostEO'};
end

statCalc = input(['Insert \n1 -> rmANOVA subsessions, \n2 -> rmANOVA segm x subsesh' ...
            '\n3 -> Pearson`s corr coeff Cz, \n4 -> FOOOF-extracted Pearson' ...
            '\n5 -> pre vs post \n6 -> pearsons at all locs \n7 -> rmANOVA Hurst' ... 
            '\n8 -> spindle correlations: ']);

switch statCalc
    case 1
        % run rmANOVA to see how power in frequency bands evolves during
        % the session
        fbands = {'theta','alpha','betaL','betaH'};
        tr_reflex = [5 6 7]; % first one is Pre
        sub_n = length(PW_all); trial_n = length(trials_nf);
        subs = repelem(1:sub_n, trial_n);
        time = repmat((1:length(trials_nf))', sub_n, 1);

        for f_ix = 1:length(fbands)
            pow_in      = zeros(length(trials_nf),sub_n);
            reflex_tr   = pow_in;
            pow_base    = bandpow.(fbands{f_ix}).PreEO(:,chix);
            for trialix = 1:length(trials_nf)
                relpow_tr = bandpow.(fbands{f_ix}).(trials_nf{trialix})(:,chix)./pow_base;
                pow_in(trialix,:) = relpow_tr;
                if ismember(trialix,tr_reflex)
                    reflex_tr(trialix,:) = 1;
                end
            end
            [p(f_ix,:), tbl, stats] = anovan(pow_in(:), {subs, time}, 'model', 'interaction', 'varnames', {'Individual', 'Time'},'continuous',2);
%             [p(f_ix,:), tbl, stats] = anovan(pow_in(:), {time, reflex_tr(:)}, 'model', 'full', 'varnames', { 'Time','Reflex'},'continuous',1);
%             [p(f_ix,:), tbl, stats] = anovan(pow_in(:), {subs, time, reflex_tr(:)}, 'model', 'full', 'varnames', {'Individual', 'Time','Reflex'},'continuous',2);
        end
    case 2
        load([folder_results 'BandPower_segm.mat']);
        % reshape struct to make table for rmANOVA
        
        % Initialize arrays to hold the variables
        numSegments = 4;
        numSubjects = length(PW_all);
        numSessions = length(trials_nf);
        Power = [];
        Session = [];
        Segment = [];
        Subject = [];
        Type    = [];
        
        % Loop over the sessions in the structure
        for i = 1:numSessions
            % Get the power values for this session
            sessionData = bandpow.alpha.(trials_nf{i})(:,:,chix);
            
%             sessionType = sessionType;
            if ismember(i,[5,6,7]) % I included baseline
                sessionType = 'r';
%                 sessionType = ones(size(sessionType));
            else
                sessionType = 'n';
%                 sessionType = zeros(size(sessionType));
            end
            % Reshape the data and append it to the Power array
            Power = [Power; sessionData(:)];
            Type = [Type; repmat(sessionType, numSegments * numSubjects, 1)];
            
            % Create arrays for the session, segment, and subject identifiers and append them to the respective variables
            Session = [Session; repmat(i, numSegments * numSubjects, 1)];
            Segment = [Segment; repmat((1:numSegments)', numSubjects, 1)];
            Subject = [Subject; reshape(repmat(1:numSubjects, numSegments, 1), [], 1)];
        end
        
        % Create a table from the arrays
        BandpowData = table(Subject, Session, Segment, Type, Power);
        
        % Write the table to a CSV file
        writetable(BandpowData, ['BandpowData_1to3_up.csv'])
        
        disp('Now it`s saved, move over to R for analysis')
        
    case 3
        % as described in RosGruzelier2011, including baseline as segm1
        
        for idxfn = 1:length(trials_nf)
            abs_alpha(:,idxfn) = bandpow.alpha.(trials_nf{idxfn})(:,chix);
            rel_alpha(:,idxfn) = abs_alpha(:,idxfn)./bandpow.alpha.PreEO(:,chix);
            
            abs_theta(:,idxfn) = bandpow.theta.(trials_nf{idxfn})(:,chix);
            abs_betaH(:,idxfn) = bandpow.betaH.(trials_nf{idxfn})(:,chix);
            
            rel_theta(:,idxfn) = bandpow.theta.(trials_nf{idxfn})(:,chix)./bandpow.theta.PreEO(:,chix);
            rel_betaH(:,idxfn) = bandpow.betaH.(trials_nf{idxfn})(:,chix)./bandpow.betaH.PreEO(:,chix);
        end
        
        training_period = 1:length(trials_nf);
% to include baseline in correlation:
%         training_period = 1:length(all_trials)-1;
%         aa = abs_alpha; abs_alpha(:,1) = aa(:,9); abs_alpha(:,2:8) = aa(:,1:7); abs_alpha(:,9) = [];
%         aa = rel_alpha; rel_alpha(:,1) = aa(:,9); rel_alpha(:,2:8) = aa(:,1:7); rel_alpha(:,9) = [];

        % Corr coeff with subsession number
        for subidx = 1:length(Subj_code)
            obs1 = rel_alpha(subidx,:);
            [pears(subidx).rho,pears(subidx).pval] = corr(obs1',training_period', 'type', 'Pearson');
            obs2 = abs_alpha(subidx,:);
            [spear(subidx).rho, spear(subidx).pval] = corr(obs2',training_period', 'type', 'Spearman');
            
            % Corr coeff beta-theta
            [thetaBetaCorrP(subidx).rho, thetaBetaCorrP(subidx).pval] = corr(rel_theta(subidx, :)', rel_betaH(subidx, :)');
            [thetaBetaCorrS(subidx).rho, thetaBetaCorrS(subidx).pval] = corr(abs_theta(subidx, :)', abs_betaH(subidx, :)', 'type', 'Spearman');
        end

        case 4
        % as described in RosGruzelier2011, including baseline as segm1
        addpath('C:\Users\2146112s\Documents\PhDina\DataAnalysis-preproc\FOOOFanalysis\')
        run step2_loadFOOOFgroup.m
                
        training_period = 1:length(trials_nf);

        % Corr coeff with subsession number
        for subidx = 1:length(Subj_code)
            obs1 = cell2mat(a_change(subidx,training_period));
            if length(obs1) == length(training_period)
                [pears(subidx).rho,pears(subidx).pval] = corr(obs1',training_period', 'type', 'Pearson');
            end
%             obs2 = abs_alpha(subidx,:);
%             [spear(subidx).rho, spear(subidx).pval] = corr(obs2',training_period', 'type', 'Spearman');
        end
        
    case 5
        abs_alpha_pre = bandpow.alpha.PreEO(:,chix);
        abs_alpha_post = bandpow.alpha.PostEO(:,chix);
        
        p_all = signrank(abs_alpha_pre,abs_alpha_post); % data not normally distributed
        % excluding unsuccessful
%         mdt = split(modType,'\');
        nfix = fct_getBestSubsession(modType,10,'abs');
        p_s = signrank(abs_alpha_pre(~isnan(nfix)),abs_alpha_post(~isnan(nfix)));
        
        % get relative change post:
        aa = 100*(abs_alpha_post./abs_alpha_pre -1);
        aa_s = 100*(abs_alpha_post(~isnan(nfix))./abs_alpha_pre(~isnan(nfix)) -1);
    case 6
        % get corr coeff at all locations on scalp 
        training_period = 1:length(trials_nf);
        
        for subidx = 1:length(Subj_code)
            for idxfn = 1:length(trials_nf)
                abs_alpha(idxfn,:) = bandpow.alpha.(trials_nf{idxfn})(subidx,:);
                rel_alpha(idxfn,:) = abs_alpha(idxfn,:)./bandpow.alpha.PreEO(subidx,:);
            end
            
            for ch_i = 1:length(abs_alpha)
                obs1 = abs_alpha(:,ch_i);
                [pears.rho(subidx,ch_i),pears.pval(subidx,ch_i)] = corr(obs1,training_period', 'type', 'Spearman');
            end
        end
    case 7
        % Hurst rmANOVA
        load(['C:\Users\2146112s\OneDrive - University of Glasgow\PhDina\healthy_study\analysis\HurstExp\H_all_' modType '.mat'])

        hurst = H_exp.H(:,1:length(trials_nf))';  % H is in the same order as trials array
        sub_n = length(PW_all); trial_n = length(trials_nf);
        subs = repelem(1:sub_n, trial_n);
        time = repmat((1:length(trials_nf))', sub_n, 1);
        reflex_tr = repmat([0 0 0 0 1 1 1 0], 1,sub_n);
        success = repelem({'U'}, 1,sub_n*trial_n);
        success(ismember(subs,success_ix)) = {'S'};
        tbl_hurst = createAnovaTable(subs, {'H-exp','Session','Reflex','Success'},hurst,time,reflex_tr,success);

        if strcmp(modType,'updown')
            n_up = find(strcmp({PW_all.group},'up'));
            group = repelem({'up'}, 1,length(n_up)*trial_n);
            group = [group repelem({'down'}, 1,(length(PW_all) - length(n_up))*trial_n)];
            tbl_hurst.group = group';
        end
        [p_h, tbl, stats] = anovan(hurst(:), {subs, time}, 'model', 'interaction', 'varnames', {'Individual', 'Time'},'continuous',2);
        [p_tr, tbl, stats] = anovan(hurst(:), {subs, time, reflex_tr(:)}, 'model', 'full', 'varnames', {'Individual', 'Time','Reflex'},'continuous',2);
        
        writetable(tbl_hurst, ['Hurst_' modType '.csv'])
        
        
        h_pre = H_exp.H(:,1);
        h_nf = mean(H_exp.H(:,2:4),2);
        h_nfall = mean(H_exp.H(:,2:8),2);
        [h,p_nf,ci,stats] = ttest(h_pre,h_nf,'alpha',0.05);
        [h,p_nfall,ci,stats] = ttest(h_pre,h_nfall,'alpha',0.05);
        
        [h,p_s,ci,stats] = ttest(H_exp.H(success_ix,1),H_exp.H(success_ix,end),'alpha',0.05);
    case 8
        % Spindle correlations
        addpath 'C:\Users\2146112s\Documents\PhDina\matlab functions, addons\swtest'
        load('C:\Users\2146112s\Documents\PhDina\healthy_data\NF learning indices\spindles_alphaup_2med.mat')
        all_trials=fieldnames(bandpow.alpha);
        all_trials(2:9) = all_trials(1:8); all_trials(1) = {'PreEO'};
        for idxfn = 1:length(all_trials)
            %     bandpow.alpha.(fn_abspow{idxfn})(idx_excl,:) = [];
            abs_alpha(idxfn,:) = bandpow.alpha.(all_trials{idxfn})(:,chix)';
            rel_alpha(idxfn,:) = bandpow_rel.alpha.(all_trials{idxfn})(:,chix)';
            % abs_theta = bandpow.theta(fixd_ord,:)';
        end
        %
        alpha_change = abs_alpha;
        alpha_change = (alpha_change./alpha_change(1,:)-1)*100; alpha_change(1,:) = [];
        alpha_allchange = mean(alpha_change(1:7,:),1);
        nf_decr = find(alpha_allchange<-10);
        
        disp('Baseline spindle IR - average alpha power NF:')
        obs1 = abs_alpha(1,:);
        obs2 = mean(NF_spin2.spin_ir(:,2:8),2);
        siglvl = 0.05;
        
        [H, pValue1, W] = swtest(obs1, siglvl);
        [H, pValue2, W] = swtest(obs2, siglvl);
        
        if pValue1 > siglvl & pValue2 > siglvl
            disp('Data is normal, you can use Pearson:')
            
            [rho,pval] = corr(obs1',obs2, 'type', 'Pearson')
        else
            disp('Spearman:')
            [rho,pval] = corr(obs1',obs2, 'type', 'Spearman')
        end
        
        figure; scatter(obs1,obs2,75,'MarkerEdgeColor','#A2142F', 'MarkerFaceColor','#A2142F')
        xlabel('\alpha power change NF (Avg)'); ylabel('PreEO NF spindle IR'); title('Baseline alpha - average change spindle IR (NF)')
        
        
        disp('Baseline spindle IR - average alpha power change NF:')
        obs1 = rel_alpha(1,:);
        obs2 = mean(NF_spin2.spin_ir(:,2:8),2);
        siglvl = 0.05;
        
        [H, pValue1, W] = swtest(obs1, siglvl);
        [H, pValue2, W] = swtest(obs2, siglvl);
        
        if pValue1 > siglvl & pValue2 > siglvl
            disp('Data is normal, you can use Pearson:')
            [rho,pval] = corr(obs1',obs2, 'type', 'Pearson')
        else
            disp('Spearman:')
            [rho,pval] = corr(obs1',obs2, 'type', 'Spearman')
        end
        
        figure; scatter(obs1,obs2,75,'MarkerEdgeColor','#A2142F', 'MarkerFaceColor','#A2142F')
        xlabel('\alpha RA PreEO'); ylabel('Change in NF spindle IR'); title('Baseline RA - spindle IR NF7')
        
end


