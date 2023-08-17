addpath('C:\Users\2146112s\Documents\PhDina\matlab functions, addons\MFDFA4\')
addpath('C:\Users\2146112s\Documents\PhDina\DataAnalysis-preproc\common_fct\')
fs = 256;
smoothwin = floor(fs/5);    % 200ms
% load data
load('C:\Users\2146112s\OneDrive - University of Glasgow\PhDina\healthy_study\data\allData_preproc\data_clean\up\session3\subHEM01_NF_7.mat')

% filter signal in frequency band of intererst & get signal envelope
[sig_env env phs] = fct_envelope_calc(EEGdata(28,:), fs, [8 12], smoothwin);

% get DFA signal
% scale defines the sample size of the non-overlapping segments in which
% the local RMS, RMS{1}, are computed
m=1;
scale = fs*[1:20];
[Ht,Htbin,Ph,Dh] = MFDFA2(sig_env,scale,m,1);

% DFA
% convert data into profile (walk)
% I think it's also mentioned in Samek2016
X=cumsum(sig_env-mean(sig_env)); 

for ns=1:length(scale)
    segments(ns)=floor(length(X)/scale(ns));
    for v=1:segments(ns)
        Idx_start=((v-1)*scale(ns))+1;
        Idx_stop=v*scale(ns);
        Index{v,ns}=Idx_start:Idx_stop;
        X_Idx=X(Index{v,ns});
        C=polyfit(Index{v,ns},X(Index{v,ns}),m);
        fit{v,ns}=polyval(C,Index{v,ns});
        RMS{ns}(v)=sqrt(mean((X_Idx-fit{v,ns}).^2));
    end
    F(ns)=sqrt(mean(RMS{ns}.^2));
end


C=polyfit(log2(scale),log2(F),1);
H1=C(1); % this is the Hurst exponent
RegLine=polyval(C,log2(scale));
figure
scatter(scale,F,'o')
hold on
refline
ylabel('log(F(time))')
xlabel('log scale (segment sample size)')