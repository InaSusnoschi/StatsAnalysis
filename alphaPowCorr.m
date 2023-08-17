% run stats corr first
a_chg = 100*(abs_alpha./abs_alpha(:,1)-1);
a_change = mean(a_chg(:,2:end),2);

b_content = bandpow.alpha.PreEO(:,28)./bandpow.broad.PreEO(:,28);
% alpha up
% a_change = [26	54	-8	18	-2	-15	-1	-31	25	-25	-40	52	-38	123	-13	-28	-1	-64	75	16];
% a_change_all = [18	54	-3	16	-3	-20	-9	-37	18	-24	-34	63	-36	111	-14	-33	-3	-64	70	7]; % incl reflex trials
% b_content= [18	18	22	67	30	25	21	20	23	55	38	10	45	21	20	41	23	54	24	19];

mdl = linregmld_plot(a_change,b_content,'change vs baseline content',1)

% theta-beta 
thetabeta = {thetaBetaCorrS.rho}';
mdl = linregmld_plot(a_change',cell2mat(thetabeta),'Alpha UP theta-beta ratio vs baseline content',1)


% alpha down
a_change = [-58	-30	-38	-59	-63	-12	16	-9	-29	-16	-15	-28	-54	-16	4	-57	22	-81	-3];
b_content= [60	34	56	40	45	66	10	24	18	21	37	33	65	34	28	45	16	65	14];

mdl = linregmld_plot(a_change',b_content','Alpha DOWN change vs baseline content',1)


% both together
a_change = [47	54	-8	18	-2	-15	-1	-31	25	-25	-40	52	-38	123	-13	-28	-1	-64	75	16 -58	-30	-38	-59	-63	-12	16	-9	-29	-16	-15	-28	-54	-16	4	-57	22	-81	-3];
b_content= [18	18	22	67	30	25	21	20	23	55	38	10	45	21	20	41	23	54	24	19 60	34	56	40	45	66	10	24	18	21	37	33	65	34	28	45	16	65	14];

mdl = linregmld_plot(a_change',b_content','Alpha BOTH change vs baseline content',1)


%% using FOOOF-extracted comps
addpath('C:\Users\2146112s\Documents\PhDina\DataAnalysis-preproc\FOOOFanalysis\')
run step2_loadFOOOFgroup.m

% only including NF 1-3
a_change_f  = [55.8 147.3 -12.9 -4.3 -38.4 -3.3 4.9 -27.7 -34.6 32.6 -19.1 -59.8 20.9 -40.6 63.2 -3.6 -38.1 41.8 46.7 -69.0 124.6 48.2];
b_content_f = [0.244320736874639;0.158827128990322;0.462697931181288;1.33723795851046;0.603193050953967;0.425465844621131;0.293648312867450;0.397037237379280;0.276010063269814;0.313295841353061;0.902964198459043;0.614205070731114;0.154016277575770;0.660336788709993;0.447809275402374;0.500349262180967;0.776729580428721;0.263193325447982;0.223799980185024;0.854325126065591;0.211850473182712;0.229889854210585];
mdl = linregmld_plot(a_change_f',a_base,'FOOOF change vs baseline power',1) 

%including reflex subsessions and NF7
a_change_all_f = [11 178	-13	-4	-46	-2	-21	-24	-40	12	-21	-39	52	-30	63	-13	-36	40	42	-64	84	18];
mdl = linregmld_plot(a_change_all_f',b_content_f,'FOOOF change vs baseline content',1) 
