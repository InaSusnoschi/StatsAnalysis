addpath('C:\Users\2146112s\Documents\PhDina\matlab functions, addons\MFDFA4\')
addpath('C:\Users\2146112s\Documents\PhDina\DataAnalysis-preproc\common_fct\')
fs = 256;
smoothwin = floor(fs/5);    % 200ms
cz_ix = 28;
% scale defines the sample size of the non-overlapping segments in which
% the local RMS, RMS{1}, are computed
% m defines polynomial order for detrending (1=linear, 2=quadratic etc.)
m=1; % NEED TO FIND OUT WHAT VALUE TO USE 
scale = fs*[1:30];

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

s_no = 3; % only calculate on session 3 data for now

folder_results  = (['C:\Users\2146112s\OneDrive - University of Glasgow\PhDina\healthy_study\analysis\spectral\' modType '\session' int2str(s_no) filesep]);
folder_data     = (['C:\Users\2146112s\OneDrive - University of Glasgow\PhDina\healthy_study\data\allData_preproc\data_clean\' modType '\session'  int2str(s_no) filesep]);
load([folder_results 'Power_allSub.mat']);

Subj_code = {PW_all.code};

% initialise matrices for H-exp, F and C
H_all = zeros(length(Subj_code),9);
F_all = zeros(length(Subj_code),9,length(scale));
C_all = zeros(length(Subj_code),9,2);

for sub_ix = 1:length(Subj_code)
    IAb = PW_all(sub_ix).IAF;
   % get files for sub
    trials = dir([folder_data Subj_code{sub_ix} '_*']);
    % move PreEO to the beginning
    tr=trials;
    trials(2:end) = tr(1:end-1); trials(1) = tr(end); clear tr
%     figure; hold on
    for tr_ix = 1:length(trials)
        load([trials(tr_ix).folder filesep trials(tr_ix).name])
        
        % filter signal in frequency band of intererst & get signal envelope
        [sig_env env phs] = fct_envelope_calc(EEGdata(cz_ix,:), fs, IAb, smoothwin);
        
        
        % DFA
        % convert data into profile (walk)
        % I think it's also mentioned in Samek2016
        X=cumsum(sig_env-mean(sig_env));
        
        [H, F, C] = getHurst(X, scale, m);
        
        % put in matrix with all sub
        H_all(sub_ix,tr_ix)     = H;
        F_all(sub_ix,tr_ix,:)   = F;
        C_all(sub_ix,tr_ix,:)   = C;
        
        RegLine=polyval(C,log2(scale));
%         scatter(scale,F,'o')
%         refline

    end
%     ylabel('log(F(time))')
%     xlabel('log scale (segment sample size)')
%     legend({'PreEO','NF1','NF2','NF3','NF4','NF5','NF6','NF7'}, 'Location','southoutside','NumColumns',9)
%     title(Subj_code{sub_ix})
end

H_exp.H = H_all;
H_exp.F = F_all;
H_exp.C = C_all;
save(['C:\Users\2146112s\OneDrive - University of Glasgow\PhDina\healthy_study\analysis\HurstExp\H_all_' modType '.mat'],'H_exp');

% plotting pre-post boxplots
figure;
boxplot([H_all(:,1),H_all(:,9)],'Labels',{'pre','post'})
hold on
allData={H_all(:,1);H_all(:,9)};
xCenter = 1:numel(allData); 
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(allData)
    plot(rand(size(allData{i}))*spread -(spread/2) + xCenter(i), allData{i}, 'mo','linewidth', 2)
end
ylabel('Hurst exponent')

% t-test in all participants
hpre = H_exp.H(:,1);
hpost = H_exp.H(:,9);
[h,p,ci,stats] = ttest(hpre,hpost);

% t-test in successful participants
nfix = fct_getBestSubsession(modType,10,'abs');
success_ix = find(~isnan(nfix));
hpre = H_exp.H(success_ix,1);
hpost = H_exp.H(success_ix,9);
[h,p_s,ci,stats] = ttest(hpre,hpost);
