load('C:\Users\2146112s\OneDrive - University of Glasgow\PhDina\iSCI\analysis\spectral\session5\relPowChange_IAF.mat')
load('C:\Users\2146112s\OneDrive - University of Glasgow\PhDina\iSCI\analysis\spectral\session5\BandPower_IAF.mat')
trialNames = fieldnames(bandpow.theta);
abspow_change = struct();
chanNo_relChange = 32;
for i = 1:numel(fieldnames(bandpow.alpha))
    abspow_change.theta(:,i)     = 100*bandpow.theta.(trialNames{i})(:,chanNo_relChange)./bandpow.theta.PreEO(:,chanNo_relChange) -100;
    abspow_change.alphaLo(:,i)   = 100*bandpow.alphaLo.(trialNames{i})(:,chanNo_relChange)./bandpow.alphaLo.PreEO(:,chanNo_relChange)  -100;
    abspow_change.alphaHi(:,i)   = 100*bandpow.alphaHi.(trialNames{i})(:,chanNo_relChange)./bandpow.alphaHi.PreEO(:,chanNo_relChange)  -100;
    abspow_change.alpha(:,i)     = 100*bandpow.alpha.(trialNames{i})(:,chanNo_relChange)./bandpow.alpha.PreEO(:,chanNo_relChange)  -100;
    abspow_change.SMR(:,i)       = 100*bandpow.SMR.(trialNames{i})(:,chanNo_relChange)./bandpow.SMR.PreEO(:,chanNo_relChange)  -100;
    abspow_change.betaH(:,i)     = 100*bandpow.betaH.(trialNames{i})(:,chanNo_relChange)./bandpow.betaH.PreEO(:,chanNo_relChange)  -100;
end


subcode = {'PCJG01'};
sub = cellstr(repmat(subcode,length(trialNames),1));
subsess = trialNames;
powchange_abs = abspow_change.alpha(2,:)';
powchange_rel = relpow_change.alpha(2,:)';

Meas = table([1:length(trialNames)]','VariableNames',{'Measurements'});

% powchange_alpha = table(subcode,powchange_abs','VariableNames',{'subject',trialNames{1:end}});
powchange_alpha = table(subcode,powchange_abs(1,:),powchange_abs(2,:),powchange_abs(3,:),powchange_abs(4,:),powchange_abs(5,:),powchange_abs(6,:),powchange_abs(7,:),'VariableNames',{'subject','meas1','meas2','meas3','meas4','meas5','meas6','meas7'});

modelspec = 'meas1-meas7~subject';

rm = fitrm(powchange_alpha,modelspec,'WithinDesign',Meas);