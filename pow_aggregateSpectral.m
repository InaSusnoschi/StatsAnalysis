addpath('C:\Users\2146112s\Documents\PhDina\matlab functions, addons\plot_topography')

folder_preproc = 'C:\Users\2146112s\Documents\PhDina\healthy_data\'; 
folder_laplace = 'C:\Users\2146112s\Documents\PhDina\healthy_data\CSDdata_onClean\';
folder_results = 'C:\Users\2146112s\OneDrive - University of Glasgow\PhDina\healthy_study\analysis\spectral\';
% channel info
load('C:\Users\2146112s\Documents\PhDina\experimental_study\z_allChanSetup\chanLocs_60.mat');
chan_labels = {chanlocs.labels};
nb_chan = length(chanlocs);
% experiment
expGroup = input('Insert 1 for alpha-UP group, 2 for alpha-DOWN group: ');
switch expGroup
    case 1
        modType = 'up';
    case 2
        modType = 'down';
    otherwise
        modType = [];
end

load([folder_results modType '\session3\Power_allSub.mat']);
load([folder_results modType '\session3\BandPower_IAF.mat']);
% participants
Subj_code = {PW_all.code}.';
trialNames = fieldnames(PW_all);
trialNames = trialNames(4:12);
trialNames(2:9) = trialNames(1:8); trialNames(1) = {'PreEO'};
fw = PW_all(1).fw; % all fw arrays are the same
% electrode grouping
af = {'AF3';'AFZ';'AF4';'F3';'F1';'FZ';'F2';'F4'};
fc = {'FC3';'FC1';'FCZ';'FC2';'FC4';'C3';'C1';'CZ';'C2';'C4'}';
% fc = {'F3';'F1';'FZ';'F2';'F4';'FC3';'FC1';'FCZ';'FC2';'FC4';'C3';'C1';'CZ';'C2';'C4'}';
cp = {'CP3';'CP1';'CPZ';'CP2';'CP4';'P3';'P1';'PZ';'P2';'P4'}';
po = {'PO3';'POZ';'PO4';'O1';'OZ';'O2'}';

fc_ix = cell2mat(cellfun(@(a) strmatch(a,chan_labels),fc,'uniform',false));
cp_ix = cell2mat(cellfun(@(a) strmatch(a,chan_labels),cp,'uniform',false));
po_ix = cell2mat(cellfun(@(a) strmatch(a,chan_labels),po,'uniform',false));
af_ix = cell2mat(cellfun(@(a) strmatch(a,chan_labels),af,'uniform',false));

PW_agg.af = zeros(length(Subj_code),length(trialNames),length(PW_all(1).fw));
PW_agg.fc = PW_agg.af;
PW_agg.cp = PW_agg.fc; PW_agg.po = PW_agg.fc;
groups = fieldnames(PW_agg);

% fixed bands
fwix_wide = fct_getFreqIx([2 30],fw);
fwix_theta = fct_getFreqIx([4 7],fw);

for subix = 1:length(Subj_code)
    IAF = PW_all(subix).IAF;
    fw_ix = fct_getFreqIx(IAF,fw);
    
    trial_ord=1;
    for trix = 1:length(trialNames)
        mat_curr = PW_all(subix).(trialNames{trix});
        if ~isempty(mat_curr)
            PW_agg.af(subix,trial_ord,:) = mean(mat_curr(:,af_ix),2);
            PW_agg.fc(subix,trial_ord,:) = mean(mat_curr(:,fc_ix),2);
            PW_agg.cp(subix,trial_ord,:) = mean(mat_curr(:,cp_ix),2);
            PW_agg.po(subix,trial_ord,:) = mean(mat_curr(:,po_ix),2);
        end
        for gr_ix = 1:length(groups)
            thetaPW.(groups{gr_ix})(subix,trial_ord) = bandpower(squeeze(PW_agg.(groups{gr_ix})(subix,trial_ord,fwix_theta(1):fwix_theta(2))),...
                fw(fwix_theta(1)+1:fwix_theta(2)+1),'psd');
            alphaPW.(groups{gr_ix})(subix,trial_ord) = bandpower(squeeze(PW_agg.(groups{gr_ix})(subix,trial_ord,fw_ix(1):fw_ix(2))),...
                fw(fw_ix(1)+1:fw_ix(2)+1),'psd');
            widePW.(groups{gr_ix})(subix,trial_ord) = bandpower(squeeze(PW_agg.(groups{gr_ix})(subix,trial_ord,fwix_wide(1):fwix_wide(2))),...
                fw(fwix_wide(1)+1:fwix_wide(2)+1),'psd');
        end
        
        trial_ord = trial_ord+1;
    end
end

best_ix_abs = zeros(length(Subj_code),length(groups));
best_ix_rel = best_ix_abs;
for gr_ix = 1:length(groups)
    alphaPW_change.(groups{gr_ix}) = 100*(alphaPW.(groups{gr_ix})(:,2:9)./alphaPW.(groups{gr_ix})(:,1)-1);
    
    alphaPW_rel.(groups{gr_ix}) = alphaPW.(groups{gr_ix})./widePW.(groups{gr_ix});
    alphaPW_rel_change.(groups{gr_ix}) = 100*(alphaPW_rel.(groups{gr_ix})(:,2:9)./alphaPW_rel.(groups{gr_ix})(:,1)-1);
    %     figure; plot(alphaPW_change.(groups{gr_ix})'); title(groups{gr_ix})
    best_ix_abs(:,gr_ix) = min(alphaPW_change.(groups{gr_ix})(:,1:7),[],2);
    best_ix_rel(:,gr_ix) = min(alphaPW_rel_change.(groups{gr_ix})(:,1:7),[],2);
end
best_ix_abs(best_ix_abs>-10) = 0;
best_ix_rel(best_ix_rel>-10) = 0;

% 
for i = 1:8
    alphaPW_change_cz(:,i)     = 100*(bandpow.alpha.(trialNames{i})(:,28)./bandpow.alpha.PreEO(:,28)  -1);
end
best_ix_cz = max(alphaPW_change_cz(:,1:7),[],2);