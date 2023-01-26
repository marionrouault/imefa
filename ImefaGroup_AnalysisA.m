% Marion Rouault IMEFA study. 2021-2023.

% This is associated with the following paper:
% Rouault, M. et al. Interoceptive and metacognitive facets of
% fatigue in multiple sclerosis. medRxiv 2023.01.23.23284429 (2023).

% This program looks at the behaviour of MS patients in relation to a
% number of psychological and physiological variables and test our
% pre-registered hypothesis A




function [] = ImefaGroup_AnalysisA()

% 1 is mfis, 2 is fss
which_fa_quest = 1;
% first component (1) or first two
% components (2) or all six regressors (no components) (6) or
% or all six regressors except for MAIA (7)
how_many_pca_compo = 1;
% colors
rgb = [[0.5,0,0.8];[0.5,0.8,0];[0.8,0.5,0]];
% display PCA subanalysis?
plot_pca = false;



% -------------- retrieve RedCap data for inter-individual analyses:

table = read_RedCap("/Users/marion/Desktop/PDOC/fatigueZurich/fa_analysis/2IMEFA-InitialVar_DATA_LABELS_2021-06-23_1109.csv");

table = table2array(table) ;

% extract data from Metacognition task 2 (ses2)
ses = double(table(:,2));

% => 71 participants have demographic data
ppidstring = table(ses==1,1);
agestring = table(ses==1,3);
sexstring = table(ses==1,4);
medistring = table(ses==1,5);
for i = 2:length(sexstring)
    ppiddemo(i) = str2num(ppidstring(i));
    age(i) = str2num(agestring(i));
    if sexstring(i)=="female"
        sex(i) = 2;
    elseif sexstring(i)=="male"
        sex(i) = 1;
    end
    if medistring(i)=="Yes"
        medi(i) = 1;
    elseif medistring(i)=="No"
        medi(i) = 0;
    end
end
% we rm PPID1 from demo who has no metacog data and is a pilot:
%  => 70 participants processed
ppiddemo = ppiddemo(2:end);
age = age(2:end);
sex = sex(2:end);
medi = medi(2:end);

% two medication regressors: whether the person takes
% immuno drugs and / or sedating drugs

medi_details_data = read_medi_details("/Users/marion/Desktop/PDOC/fatigueZurich/fa_analysis/Imefa_medication_ZMM_220717.csv");

medi_immuno = medi_details_data(ses==1,3);
medi_immuno = table2array(medi_immuno);
medi_immuno = medi_immuno(2:end);

medi_sedate = medi_details_data(ses==1,9);
medi_sedate = table2array(medi_sedate);
medi_sedate = medi_sedate(2:end);



% -------------- retrieve RedCap data for sleep metrics

sleep_data = read_sleep("/Users/marion/Desktop/PDOC/fatigueZurich/fa_analysis/2IMEFA-Sleep_DATA_2022-04-04_1511.csv");

sleep_data_labels = read_sleep("/Users/marion/Desktop/PDOC/fatigueZurich/fa_analysis/2IMEFA-Sleep_DATA_LABELS_2022-04-04_1511.csv");

% PSQI Sleep quality index questionnaire, first question is column13
% psqi total score is the sum of seven components:

ppidpsqi = table2array(sleep_data_labels(ses==2,1));

% Component 1: Subjective sleep quality?question 9
psqi_cmp1 = sleep_data(ses==2,31);
psqi_cmp1 = table2array(psqi_cmp1)';

% Component 2: Sleep latency?questions 2 and 5a
psqi_q2 = sleep_data(ses==2,14);

psqi_q2 = double(table2array(psqi_q2));

% convert score of Q2
for j = 1:length(psqi_q2)
    if psqi_q2(j) <= 15
        psqi_q2(j) = 0;
    elseif (psqi_q2(j) > 15) && (psqi_q2(j) <= 30)
        psqi_q2(j) = 1;
    elseif (psqi_q2(j) > 30) && (psqi_q2(j) <= 60)
        psqi_q2(j) = 2;
    elseif psqi_q2(j) > 60
        psqi_q2(j) = 3;
    end
end
clear j

psqi_q5a = sleep_data(ses==2,17);
psqi_q5a = table2array(psqi_q5a);

for j = 1:length(psqi_q2)
    if psqi_q2(j) + psqi_q5a(j) == 0
        psqi_cmp2(j) = 0;
    elseif (psqi_q2(j) + psqi_q5a(j) == 1) || (psqi_q2(j) + psqi_q5a(j) == 2)
        psqi_cmp2(j) = 1;
    elseif (psqi_q2(j) + psqi_q5a(j) == 3) || (psqi_q2(j) + psqi_q5a(j) == 4)
        psqi_cmp2(j) = 2;
    elseif psqi_q2(j) + psqi_q5a(j) > 4
        psqi_cmp2(j) = 3;
    end
end
clear j

% Component 3: Sleep duration?question 4
psqi_q4 = sleep_data(ses==2,16);

psqi_q4 = double(table2array(psqi_q4));

for j = 1:length(psqi_q4)
    if psqi_q4(j) >= 7
        psqi_cmp3(j) = 0;
    elseif (psqi_q4(j) >= 6) && (psqi_q4(j) < 7)
        psqi_cmp3(j) = 1;
    elseif (psqi_q4(j) >= 5) && (psqi_q4(j) < 6)
        psqi_cmp3(j) = 2;
    elseif psqi_q4(j) < 5
        psqi_cmp3(j) = 3;
    end
end
clear j


% Component 4: Sleep efficiency?questions 1, 3, and 4
psqi_q1 = table2array(sleep_data(ses==2,13));%heure au lit
psqi_q3 = table2array(sleep_data(ses==2,15));%heure debout

% manual entry for sleep hours extraction:
hours_in_bed = [8.83 6.5 7.5 9 8 6.5 9.75 8.5 11 9.5 8.5 6.33 8.66 7.5 ...
    9 7.2 9.5 8.5 4.75 10 7.083 9 11 8.75 7.33 8 7.83 8 8.25 8.33 8.5 9.5 ...
    6.5 8 9.5 7.15 10 7.33 8 6.5 9 9.5 8.75 6.5 10 8.5 8.33 8 7.75 10 8.5 ...
    9.75 6.5 9.75 12.7 8.5 10 10 10 9.5 8 7 9.5 8.5 8.42 7.66 8 11]';

sleep_efficiency = (psqi_q4./hours_in_bed)*100;

for j = 1:length(psqi_q1)
    if sleep_efficiency(j) >= 85
        psqi_cmp4(j) = 0;
    elseif (sleep_efficiency(j) > 75) && (sleep_efficiency(j) <= 85)
        psqi_cmp4(j) = 1;
    elseif (sleep_efficiency(j) > 65) && (sleep_efficiency(j) <= 75)
        psqi_cmp4(j) = 2;
    elseif sleep_efficiency(j) < 65
        psqi_cmp4(j) = 3;
    end
end
clear j



% Component 5: Sleep disturbance?questions 5b-5j
psqi_q5bj = sum(table2array(sleep_data(ses==2,18:26)),2);
for j = 1:length(psqi_q5bj)
    if psqi_q5bj(j) == 0
        psqi_cmp5(j) = 0;
    elseif (psqi_q5bj(j) >= 1) && (psqi_q5bj(j) <= 9)
        psqi_cmp5(j) = 1;
    elseif (psqi_q5bj(j) >= 10) && (psqi_q5bj(j) <= 18)
        psqi_cmp5(j) = 2;
    elseif psqi_q5bj(j) >= 19
        psqi_cmp5(j) = 3;
    end
end
clear j


% Component 6: Use of sleep medication?question 6
psqi_q6 = table2array(sleep_data(ses==2,28));
psqi_cmp6 = psqi_q6';


% Component 7: Daytime dysfunction?questions 7 and 8
psqi_q7 = table2array(sleep_data(ses==2,29));
psqi_q8 = table2array(sleep_data(ses==2,30));

for j = 1:length(psqi_q7)
    if psqi_q7(j) + psqi_q8(j) == 0
        psqi_cmp7(j) = 0;
    elseif (psqi_q7(j) + psqi_q8(j) == 1) || (psqi_q7(j) + psqi_q8(j) == 2)
        psqi_cmp7(j) = 1;
    elseif (psqi_q7(j) + psqi_q8(j) == 3) || (psqi_q7(j) + psqi_q8(j) == 4)
        psqi_cmp7(j) = 2;
    elseif psqi_q7(j) + psqi_q8(j) > 4
        psqi_cmp7(j) = 3;
    end
end
clear j


psqi_tot = psqi_cmp1 + psqi_cmp2 + psqi_cmp3 + psqi_cmp4 + ...
    psqi_cmp5 + psqi_cmp6 + psqi_cmp7;

% filter out patients who did do psqi and align with PPID
% => 68 participants have psqi data
psqi_tot = psqi_tot(~isnan(psqi_tot))';


% -------------- retrieve RedCap data for diagnostic dates

dates = read_diseaseDUR("/Users/marion/Desktop/PDOC/fatigueZurich/fa_analysis/2IMEFA-Diseaseduration_DATA_LABELS_2021-08-03_1802.csv");

dates = dates(ses==1,:);

dates = dates(2:end,:);

consent_day = table2array(dates(:,5));
diagnostic_day = table2array(dates(:,6));

diseaseDuration = hours(consent_day-diagnostic_day);
diseaseDuration(diseaseDuration<0)=0; %one person has mistake in dates


% -------------- retrieve RedCap data for questionnaires FSS GSES MAIA

fssmaiagses = read_FFSMAIAGSES("/Users/marion/Desktop/PDOC/fatigueZurich/fa_analysis/2IMEFA-FFSMAIAGSES_DATA_2021-07-05_1812.csv");

fssmaiagses2 = fssmaiagses(ses==2,:);
for i = 1:size(fssmaiagses2,1)
    % General self-efficacy scale
    gses_tot(i) = sum(table2array(fssmaiagses2(i,14:23)));
    
    % MAIA-3 subscale Not Worrying: ((5 - Q8) + (5 - Q9) + Q10) / 3
    maia_not_worrying(i) = ((5-table2array(fssmaiagses2(i,31))) + (5-table2array(fssmaiagses2(i,32))) + table2array(fssmaiagses2(i,33)))/3 ;
    % MAIA-8 subscale Trusting: (Q30 + Q31 + Q32) / 3
    maia_trusting(i) = (table2array(fssmaiagses2(i,53)) + table2array(fssmaiagses2(i,54)) + table2array(fssmaiagses2(i,55)))/3 ;
    % Sum of subscales 3 and 8:
    maia38(i) = maia_not_worrying(i) + maia_trusting(i);
    
    % filter out patients who did do gses
    ppidgsesmaia(i) = table2array(fssmaiagses2(i,1)); % 68 participants have gses and maia data
end

fssmaiagses1 = fssmaiagses(ses==1,:);
for i = 1:size(fssmaiagses1,1)
    % Fatigue scale FSS
    fss_tot(i) = sum(table2array(fssmaiagses1(i,5:13)))/9;
    % filter out patients who did do FSS
    ppidfss(i) = table2array(fssmaiagses1(i,1)); % 71 participants have fss data
end



% MFIS has 21 items: The total MFIS score can range from 0 to 84
% MFIS was collected during session 2

tblraw = read_RedCap_raw("/Users/marion/Desktop/PDOC/fatigueZurich/fa_analysis/2IMEFA-InitialVar_DATA_2021-06-23_1352.csv");

tblraw = table2array(tblraw) ;

for i = 1:size(tblraw,1)
    mfis_tot(i) = sum(tblraw(i,2:22));
end

% filter out patients who did do mfis and align with PPID
ppidmfis = tblraw(~isnan(mfis_tot),1); 
% => 61 participants have mfis data
mfis_tot = mfis_tot(~isnan(mfis_tot));




% -------------- retrieve RedCap data for physiological measurements

InterocepvarDATA = read_homeo("/Users/marion/Desktop/PDOC/fatigueZurich/fa_analysis/2IMEFA-Interocepvar_DATA_2021-07-12_1430.csv");

ppidses2 = table2array(InterocepvarDATA(ses==2,1));

% sudomotor activity averaging sudomotor_hands and sudomotor_feet
sudo = table2array(InterocepvarDATA(ses==2,5:6));
ppidsudo = [];
for subj = 1:size(sudo,1)
    if ~isnan(sudo(subj,1)) && ~isnan(sudo(subj,2))
        ppidsudo = [ppidsudo ppidses2(subj)];
        % => 65 participants have sudomotor data
    end
end
sudo = sudo(~isnan(sudo(:,1)),:);

% HRV (heart rate variability)
% From a physiological point of view, deep breathing is most appropriate
% The most reliable indicator of parasympathetic activity would
% probably be the RMSSD (for deep breathing)
hrv = table2array(InterocepvarDATA(ses==2,11));
ppidhrv = [];
for subj = 1:length(hrv)
    if ~isnan(hrv(subj,1))
        ppidhrv = [ppidhrv ppidses2(subj)];
        % => 68 participants have hrv data
    end
end

% We are interested in autonomous regulation, therefore the most appropriate approach
% is the following difference: value immediate standing minus 10minute lying
% separately for blood pressure systolic, diastolic, and heart rate

%deltaBP lying vs. standing (blood pressure, ast_bpxxx)
%deltaHR lying vs. standing (heart rate, ast_hrxxx)
systo_10min = table2array(InterocepvarDATA(ses==2,18));
systo_imm = table2array(InterocepvarDATA(ses==2,19));

diasto_10min = table2array(InterocepvarDATA(ses==2,20));
diasto_imm = table2array(InterocepvarDATA(ses==2,21));

hr_10min = table2array(InterocepvarDATA(ses==2,22));
hr_imm = table2array(InterocepvarDATA(ses==2,23));

% Difference immediate standing minus 10min lying so that highest
% difference corresponds to best autonomic regulation:
deltaBP_systo  = systo_imm - systo_10min;
deltaBP_diasto = diasto_imm - diasto_10min;
deltaHR = hr_imm - hr_10min;

ppidheart = [];
for subj = 1:size(deltaBP_systo,1)
    if ~isnan(deltaBP_systo(subj)) && ~isnan(deltaBP_diasto(subj)) && ~isnan(deltaHR(subj))
        ppidheart = [ppidheart ppidses2(subj)];
        % => 65 participants have sudomotor data
    end
end
deltaBP_systo = deltaBP_systo(~isnan(deltaBP_systo));
deltaBP_diasto = deltaBP_diasto(~isnan(deltaBP_diasto));
deltaHR = deltaHR(~isnan(deltaHR));


% restrict to participants who have ALL homeo measures: 63 participants
ppidhomeo = intersect(ppidsudo,intersect(ppidhrv,ppidheart));



x1 = [];
x2 = [];
x3 = [];
x4 = [];
x5 = [];
for su = 1:71
    if ismember(su,ppidhomeo)
        x1 = [x1 ; deltaBP_systo(ppidheart==su)];
        x2 = [x2 ; deltaBP_diasto(ppidheart==su)];
        x3 = [x3 ; deltaHR(ppidheart==su)];
        x4 = [x4 ; sudo(ppidsudo==su,:)];
        x5 = [x5 ; hrv(ppidhrv==su,:)];
    end
end


% aparte: correlation matrix between measurements
R1 = corrcov(cov([x1 x2 x3 x4 x5]));
if plot_pca
    hf = figure('Color','white');
    heatmap(R1,'Colormap',parula,'FontSize',15,...
        'ColorLimits',[-1,1], ...
        'XData',{'deltaBP systolic','deltaBP diastolic',...
        'deltaHR','sudoHands','sudoFeet','HRV'}, ...
        'YData',{'deltaBP systolic','deltaBP diastolic',...
        'deltaHR','sudoHands','sudoFeet','HRV'});
    set(gca,'FontName','Helvetica','FontSize',15);
    figure(hf);
    fname = './FigurePCAphysio';
    print(fname,'-painters','-dpdf');
end



all_homeo = [x1 x2 x3 x4 x5];


[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(all_homeo);
% initial preregistered analysis with first component:
coef_homeo = SCORE(:,1);
% alternative analysis using first-two components:
homeo_pca1 = SCORE(:,1);
homeo_pca2 = SCORE(:,2);


if plot_pca
    disp 'PCA homeostasis. Percentage of variance explained by each component:'
    for component = 1:length(EXPLAINED)
        disp(['Component ',num2str(component),': ', ...
            num2str(EXPLAINED(component)),'%'])
    end
    
    figure;
    bar(LATENT);
    set(gca,'fontsize',20)
    ylabel('eigenvalues','fontsize',20)
    title('PCA over 6 physiological measurements','fontsize',20)
    
    figure;
    bar(COEFF(:,1:2));
    ylabel('loadings','fontsize',10)
    xlabel('physiological measurements','fontsize',10)
    set(gca,'XTick',1:6,'XTickLabel',{'deltaBP systolic','deltaBP diastolic',...
        'deltaHR','sudoHands','sudoFeet','HRV deep breathing'},'fontsize',10)
    title('PCA over 6 physio measures for 2 main components','fontsize',10)
end





% -------------- Test our pre-registered hypotheses --------------

% A. Perceived fatigue is related to measures of interoception and autonomic regulation

% - Dependent variable: FSS or MFIS scores
% - Regressors of interest:
% -- RegI1: sum of the MAIA3 [Not-Worrying] and MAIA8 [Trusting]
% -- RegI2: first principal component of autonomic measures (HRV, deltaBP and
% deltaHR (sitting vs. standing), sudomotor activity, and sympathetic skin response)
% - Regressors of no interest: duration of disease, age, sex, medication (dose and type).

% Hypothesis A1: (RegI1)
% Feeling of being in homeostasis and control (i.e. MAIA3+8) are negatively associated with FSS
% Hypothesis A2: (RegI2)
% Integrity of homeostatic regulation and autonomic function are negatively associated with FSS


% find restricted pool of participants who have demo + auton + mfis data

% initialise dependent variable (either fatigue score)
y_mfis = [];
y_fss = [];
% initialise regressors of interest
x_maia38 = []; % RegI1
x_homeo  = []; % RegI2
% alternative analysis where x_homeo is split into first-two components:
x_homeo1  = [];
x_homeo2  = [];
% alternative analysis where x_homeo is actually the 6 regressors:
x_homeoR1  = [];
x_homeoR2  = [];
x_homeoR3  = [];
x_homeoR4  = [];
x_homeoR5  = [];
x_homeoR6  = [];

% initialise regressors of no interest
x_age = [];
x_sex = [];
x_medi = [];
x_medi_immuno = [];
x_medi_sedate = [];
x_dur = [];
x_sleep = [];

if which_fa_quest == 1
    ppidfa = ppidmfis;
elseif which_fa_quest == 2
    ppidfa = ppidfss;
end

for s = 1:length(ppidstring) % the largest sample
    
    if ismember(s,ppidfa) && ismember(s,ppidgsesmaia) && ismember(s,ppidhomeo) && ismember(s,ppiddemo) && ismember(s,ppidpsqi)
        % patient s has all required types of data: include in analysis
        
        % define dependent variable (either of the two fatigue scores)
        if which_fa_quest == 1
            y_mfis = [y_mfis ; mfis_tot(ppidmfis==s)];
        elseif which_fa_quest == 2
            y_fss  = [y_fss  ; fss_tot(ppidfss==s)];
        end
        % define regressors of interest
        x_maia38 = [x_maia38  ; maia38(ppidgsesmaia==s)];
        x_homeo  = [x_homeo  ; coef_homeo(ppidhomeo==s)];
        x_homeo1 = [x_homeo1  ; homeo_pca1(ppidhomeo==s)];
        x_homeo2 = [x_homeo2  ;homeo_pca2(ppidhomeo==s)];
        x_homeoR1 = [x_homeoR1  ; x1(ppidhomeo==s)];
        x_homeoR2 = [x_homeoR2  ; x2(ppidhomeo==s)];
        x_homeoR3 = [x_homeoR3  ; x3(ppidhomeo==s)];
        x_homeoR4 = [x_homeoR4  ; x4(ppidhomeo==s,1)];
        x_homeoR5 = [x_homeoR5  ; x4(ppidhomeo==s,2)];
        x_homeoR6 = [x_homeoR6  ; x5(ppidhomeo==s)];
        
        % define regressors of no interest
        x_age  = [x_age ; age(ppiddemo==s)];
        x_sex  = [x_sex ; sex(ppiddemo==s)];
        x_medi = [x_medi ; medi(ppiddemo==s)];
        x_medi_immuno = [x_medi_immuno ; medi_immuno(ppiddemo==s)];
        x_medi_sedate = [x_medi_sedate ; medi_sedate(ppiddemo==s)];
        x_dur  = [x_dur ; diseaseDuration(ppiddemo==s)];
        x_sleep= [x_sleep ; psqi_tot(ppidpsqi==s)];
        
    else
        disp(['missing some of the data for patient #',num2str(s)])
    end
end
% => 53 patients remained in this restricted pool (MFIS) (63 for FSS)


% side note: visualise correlation matrix between measurements
R2 = corrcov(cov([y_mfis x_maia38 x_homeo x_age x_sex x_medi x_dur x_sleep]));
R3 = corrcov(cov([y_mfis x_maia38 x_homeoR1 x_homeoR2 x_homeoR3 x_homeoR4 x_homeoR5 x_homeoR6]));
R4 = corrcov(cov([y_mfis x_maia38 x_homeo x_age x_sex x_medi_immuno x_medi_sedate x_dur x_sleep]));
if plot_pca
    figure;
    heatmap(R2,'Colormap',parula,'FontSize',15,...
        'ColorLimits',[-1,1], ...
        'XData',{'MFIS','MAIA','PC1','age','sex',...
        'medic','diseaseDur','sleep'}, ...
        'YData',{'MFIS','MAIA','PC1','age','sex',...
        'medic','diseaseDur','sleep'});
    set(gca,'FontName','Helvetica','FontSize',15);
    
    figure;
    heatmap(R3,'Colormap',parula,'FontSize',15,...
        'ColorLimits',[-1,1], ...
        'XData',{'MFIS','MAIA','deltaBP systolic','deltaBP diastolic',...
        'deltaHR','sudoHands','sudoFeet','HRV deep breathing'}, ...
        'YData',{'MFIS','MAIA','deltaBP systolic','deltaBP diastolic',...
        'deltaHR','sudoHands','sudoFeet','HRV deep breathing'});
    set(gca,'FontName','Helvetica','FontSize',15);
    
    figure;
    heatmap(R4,'Colormap',parula,'FontSize',15,...
        'ColorLimits',[-1,1], ...
        'XData',{'MFIS','MAIA','PC1','age','sex',...
        'mediImmuno','mediSedate','diseaseDur','sleep'}, ...
        'YData',{'MFIS','MAIA','PC1','age','sex',...
        'mediImmuno','mediSedate','diseaseDur','sleep'});
    set(gca,'FontName','Helvetica','FontSize',15);
end


% zscore our regressors for comparability of regressions coefs:
x_maia38 = (x_maia38-mean(x_maia38))/std(x_maia38);
x_homeo = (x_homeo-mean(x_homeo))/std(x_homeo);
x_homeo1 = (x_homeo1-mean(x_homeo1))/std(x_homeo1);
x_homeo2 = (x_homeo2-mean(x_homeo2))/std(x_homeo2);
x_homeoR1 = (x_homeoR1-mean(x_homeoR1))/std(x_homeoR1);
x_homeoR2 = (x_homeoR2-mean(x_homeoR2))/std(x_homeoR2);
x_homeoR3 = (x_homeoR3-mean(x_homeoR3))/std(x_homeoR3);
x_homeoR4 = (x_homeoR4-mean(x_homeoR4))/std(x_homeoR4);
x_homeoR5 = (x_homeoR5-mean(x_homeoR5))/std(x_homeoR5);
x_homeoR6 = (x_homeoR6-mean(x_homeoR6))/std(x_homeoR6);
x_age = (x_age-mean(x_age))/std(x_age);
x_sex = (x_sex-mean(x_sex))/std(x_sex);
x_medi_immuno = (x_medi_immuno-mean(x_medi_immuno))/std(x_medi_immuno);
x_medi_sedate = (x_medi_sedate-mean(x_medi_sedate))/std(x_medi_sedate);
x_dur = (x_dur-mean(x_dur))/std(x_dur);
x_sleep = (x_sleep-mean(x_sleep))/std(x_sleep);

% build regressor matrix
if how_many_pca_compo == 2
    X = [x_maia38 x_homeo1 x_homeo2 x_age x_sex x_medi_immuno x_medi_sedate x_dur x_sleep];
elseif how_many_pca_compo == 1
    X = [x_maia38 x_homeo x_age x_sex x_medi_immuno x_medi_sedate x_dur x_sleep];
elseif how_many_pca_compo == 6 % no PCs, raw 6 reg of origin
    X = [x_maia38 x_homeoR1 x_homeoR2 x_homeoR3 x_homeoR4 x_homeoR5 x_homeoR6 ...
        x_age x_sex x_medi_immuno x_medi_sedate x_dur x_sleep];
    Xalt = [x_homeoR1 x_homeoR2 x_homeoR3 x_homeoR4 x_homeoR5 x_homeoR6];
elseif how_many_pca_compo == 7 % no PCs, raw 6 reg of origin minus MAIA
    X = [x_homeoR1 x_homeoR2 x_homeoR3 x_homeoR4 x_homeoR5 x_homeoR6 ...
        x_age x_sex x_medi_immuno x_medi_sedate x_dur x_sleep];
end

if how_many_pca_compo == 1
    disp 'Model with one PCA component for homeostasis:'
    if which_fa_quest == 1
        % perform regression
        [betaA,~,statsA] = glmfit(X,y_mfis);
        % perform GLM to get F-test
        mdl = fitglm(X,y_mfis,'y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8','Distribution','normal');
        p_mdl = coefTest(mdl);
        disp 'Ftest Model A with MFIS, pvalue='
        p_mdl
    elseif which_fa_quest == 2
        % perform regression
        [betaA,~,statsA] = glmfit(X,y_fss);
        % perform GLM to get F-test
        mdl = fitglm(X,y_fss,'y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8','Distribution','normal');
        p_mdl = coefTest(mdl);
        disp 'Ftest Model A with FSS, pvalue='
        p_mdl
    end
elseif how_many_pca_compo == 2
    disp 'Model with two PCA component for homeostasis:'
    if which_fa_quest == 1
        % perform regression
        [betaA,~,statsA] = glmfit(X,y_mfis);
        % perform GLM to get F-test
        mdl = fitglm(X,y_mfis,'y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9','Distribution','normal');
        p_mdl = coefTest(mdl);
        disp 'Ftest Model A with MFIS, pvalue='
        p_mdl
    elseif which_fa_quest == 2
        % perform regression
        [betaA,~,statsA] = glmfit(X,y_fss);
        % perform GLM to get F-test
        mdl = fitglm(X,y_fss,'y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9','Distribution','normal');
        p_mdl = coefTest(mdl);
        disp 'Ftest Model A with FSS, pvalue='
        p_mdl
    end
elseif how_many_pca_compo == 6
    disp 'Model with all six physio reg (no PCs) for homeostasis:'
    if which_fa_quest == 1
        % perform regression
        [betaA,~,statsA] = glmfit(X,y_mfis);
        % perform GLM to get F-test
        mdl = fitglm(X,y_mfis,'y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 +x10 + x11 + x12 + x13','Distribution','normal');
        p_mdl = coefTest(mdl);
        disp 'Ftest Model A with MFIS, pvalue='
        p_mdl
        p_mdl = coefTest(mdl,[0 0 1 1 1 1 1 1 0 0 0 0 0 0]);
        disp 'Ftest Model A with MFIS, test only the six physio regz, pvalue='
        p_mdl
        
        % perform regression
        [betaA,~,statsA] = glmfit(Xalt,y_mfis);
        mdl = fitglm(Xalt,y_mfis,'y ~ x1 + x2 + x3 + x4 + x5 + x6','Distribution','normal');
        p_mdl = coefTest(mdl);
        disp 'Ftest Model A with MFIS and only six physio predictors, pvalue='
        p_mdl
    elseif which_fa_quest == 2
        % perform regression
        [betaA,~,statsA] = glmfit(X,y_fss);
        % perform GLM to get F-test
        mdl = fitglm(X,y_fss,'y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 +x10 + x11 + x12 + x13','Distribution','normal');
        p_mdl = coefTest(mdl);
        disp 'Ftest Model A with FSS, pvalue='
        p_mdl
        p_mdl = coefTest(mdl,[0 0 1 1 1 1 1 1 0 0 0 0 0 0]);
        disp 'Ftest Model A with FSS, test only the six physio regz, pvalue='
        p_mdl
        
        % perform regression
        [betaA,~,statsA] = glmfit(Xalt,y_fss);
        mdl = fitglm(Xalt,y_fss,'y ~ x1 + x2 + x3 + x4 + x5 + x6','Distribution','normal');
        p_mdl = coefTest(mdl);
        disp 'Ftest Model A with FSS and only six physio predictors, pvalue='
        p_mdl
    end
elseif how_many_pca_compo == 7
    disp 'Model with all six physio reg (no PCs) minus MAIA:'
    if which_fa_quest == 1
        % perform regression
        [betaA,~,statsA] = glmfit(X,y_mfis);
        % perform GLM to get F-test
        mdl = fitglm(X,y_mfis,'y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 +x10 + x11 + x12','Distribution','normal');
        p_mdl = coefTest(mdl);
        disp 'Ftest Model A with MFIS, pvalue='
        p_mdl
        p_mdl = coefTest(mdl,[0 1 1 1 1 1 1 0 0 0 0 0 0]);
        disp 'Ftest Model A with MFIS, test only the six physio regz, pvalue='
        p_mdl
        
    elseif which_fa_quest == 2
        % perform regression
        [betaA,~,statsA] = glmfit(X,y_fss);
        % perform GLM to get F-test
        mdl = fitglm(X,y_fss,'y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 +x10 + x11 + x12','Distribution','normal');
        p_mdl = coefTest(mdl);
        disp 'Ftest Model A with FSS, pvalue='
        p_mdl
        p_mdl = coefTest(mdl,[0 1 1 1 1 1 1 0 0 0 0 0 0]);
        disp 'Ftest Model A with FSS, test only the six physio regz, pvalue='
        p_mdl
        
    end
end

disp 'Statistics for each regressor:'
disp 'betas:'
statsA.beta(2:end)
disp 'tvalues:'
statsA.t(2:end)
disp 'pvalues:'
statsA.p(2:end)


% examine Benjamini-Hochler FDR correction:

orig_pvalues = [[1:length(statsA.p(2:end))]' statsA.p(2:end)];
% order the pvalues:
ordered_pvalues = sortrows(orig_pvalues,2);
% add rank:
ordered_pvalues = [ordered_pvalues [1:length(statsA.p(2:end))]'];
% what percentage of FDR rate you decide
fdr_rate = .05;
% how many tests do you want to correct for
m = 3;
% calculate each individual p-value?s Benjamini-Hochberg critical value
ordered_pvalues(:,4) = ordered_pvalues(:,3)./m*fdr_rate;
% identify the highest p-value that is also smaller than the critical value
ordered_pvalues



% -------------- Figure Hypothesis A --------------

hf = figure('Color','white');

hold on
if which_fa_quest == 1
    ylim([-10,11]);%scale for MFIS
elseif which_fa_quest == 2
    ylim([-0.8,1.2]);%scale for FSS
end
nbar = length(betaA)-1;
pbar = 1/2*(nbar+0.2)/2.2*2;
for i = 1:nbar
    x = betaA(i+1);
    bar(i,mean(x),0.8,'EdgeColor','none','FaceColor',0.5*(rgb(1,:)+1),'FaceAlpha',0.5);
    pos = i;
    bar(pos,mean(x),0.8,'EdgeColor',rgb(1,:),'FaceColor','none','LineWidth',1);
    plot(pos*[1,1],[mean(x)-(statsA.se(i+1)) mean(x)+(statsA.se(i+1))],'k','LineWidth',1);
end
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',8);
ylabel('Regression coefficients','FontName','Helvetica','FontSize',8);
if how_many_pca_compo == 2
    set(gca,'XTick',1:size(X,2),'XTickLabel',{'MAIA_3_8','homeo_1','homeo_2','age','sex','medic','duration','sleep'});
    fname = './FigureA_PC1andPC2';
elseif how_many_pca_compo == 1
    set(gca,'XTick',1:size(X,2),'XTickLabel',{'MAIA_3_8','homeo','age','sex','mediImmuno','medicSedate','duration','sleep'});
    fname = './FigureA_PC1';
end
if which_fa_quest == 1
    xlabel('Contributors to fatigue scores (MFIS)','FontName','Helvetica','FontSize',8);
    fscale = '_MFIS';
elseif which_fa_quest == 2
    xlabel('Contributors to fatigue scores (FSS)','FontName','Helvetica','FontSize',8);
    fscale = '_FSS';
end

set(hf,'PaperPositionMode','manual', ...
    'PaperPosition',[2.5,13,16,4],'PaperUnits','centimeters', ...
    'PaperType','A4','PaperOrientation','portrait');
figure(hf);
print([fname,fscale],'-painters','-dpdf');

end
