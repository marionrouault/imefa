% Marion Rouault IMEFA study. 2021-2023.

% This is associated with the following paper:
% Rouault, M. et al. Interoceptive and metacognitive facets of
% fatigue in multiple sclerosis. medRxiv 2023.01.23.23284429 (2023).

% This program looks at the behaviour of MS patients in relation to a
% number of psychological and physiological variables and test our
% pre-registered hypothesis C: Measures of interoception and autonomic
% regulation are related to measures of metacognition



function [] = ImefaGroup_AnalysisC()

% test hypothesis C1 or C2
which_hyp = 2;
% colors
rgb = [[0.5,0.8,0];[0.8,0.5,0];[0,0.2,0.8]; ...
    [0.8,0,0.2];[0.5,0,0.8]];








% -------------- retrieve metacognition data for inter-individual analyses:

load groupSES2

ppidmeta2    = groupSES2(:,1);
acc2         = groupSES2(:,2);
metacogbias2 = groupSES2(:,4);
metacogeff   = groupSES2(:,11);



% -------------- retrieve RedCap data for inter-individual analyses:

table = read_RedCap("/Users/marion/Desktop/PDOC/fatigueZurich/fa_analysis/2IMEFA-InitialVar_DATA_LABELS_2021-06-23_1109.csv");

table = table2array(table) ;

% extract data from Metacognition task 2 (ses2)
ses = double(table(:,2));

% 71 participants have demographic data
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
diseaseDuration(diseaseDuration<0)=0; % one person has mistake in dates



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
    ppidfss(i) = table2array(fssmaiagses1(i,1)); 
    % => 71 participants have fss data
end


% -------------- retrieve RedCap data for questionnaires MSES

msse = read_MSSE("/Users/marion/Desktop/PDOC/fatigueZurich/fa_analysis/2IMEFA-MSSE_DATA_2022-04-29_1255.csv");

msse_ses2 = msse(ses==2,:);
for i = 1:size(fssmaiagses2,1)
    
    % MS self-efficacy scale
    msse_tot(i) = sum(table2array(msse_ses2(i,5:14)));
    
    % filter out patients who did do gses
    ppidmsse(i) = table2array(msse_ses2(i,1)); 
    % => 68 participants have gses data
end




% -------------- retrieve RedCap data for physiological measurements

if which_hyp == 2
    
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
    
    
    % We are interested in autonomous regulation, therefore the most appropriate approach would
    % be to take the following difference: immediate standing minus 10minute lying
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
            % => 65 participants have heart data
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
    
    all_homeo = [x1 x2 x3 x4 x5];
    
    % run the PCA over all physiological measurements:
    [COEFF, SCORE, LATENT]  = pca(all_homeo);
    
    coef_homeo = SCORE(:,1); % first component will be our regressor
    
    % plot the spectrum of the eigenvalues when doing the PCA
    % => to examine whether dimensionality reduction via PCA
    % actually makes sense for the physiological data:
    
    hf = figure('Color','white');
    hold on
    bar(LATENT);
    ylabel('eigenvalues','fontsize',15)
    title('PCA over physiological measures','fontsize',15)
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[1,1,1]);
    set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(1,1));
    set(gca,'FontName','Helvetica','FontSize',15);
    figure(hf);
    fname = './FigureImefa_PCAphysio';
    print(fname,'-painters','-dpdf');
    
end





% find restricted pool of participants who have metacog + demo + maia or homeo data
x_metacogbias2 = [];
x_metacogeff = [];
x_acc2 = [];
if which_hyp == 1
    y_maia38 = [];
elseif which_hyp == 2
    y_homeo = [];
end
x_age = [];
x_sex = [];
x_medi = [];
x_medi_immuno = [];
x_medi_sedate = [];
x_dur = [];
x_sleep = [];

x_gses = [];
x_msse = [];

for s = 1:length(ppiddemo) % the largest sample
    if which_hyp == 1
        if ismember(s,ppidgsesmaia) && ismember(s,ppidmeta2) && ismember(s,ppiddemo) && ismember(s,ppidpsqi)
            % patient s has all required types of data: include in analysis
            
            % define dependent variable
            y_maia38 = [y_maia38 ; maia38(ppidgsesmaia==s)];
            
            
            % define regressors of interest
            x_metacogbias2 = [x_metacogbias2 ; metacogbias2(ppidmeta2==s)];
            x_metacogeff = [x_metacogeff ; metacogeff(ppidmeta2==s)];
            
            % define regressors of no interest
            x_acc2 = [x_acc2 ; acc2(ppidmeta2==s)];
            
            x_age  = [x_age ; age(ppiddemo==s)];
            x_sex  = [x_sex ; sex(ppiddemo==s)];
            x_medi = [x_medi; medi(ppiddemo==s)];
            x_medi_immuno = [x_medi_immuno; medi_immuno(ppiddemo==s)];
            x_medi_sedate = [x_medi_sedate; medi_sedate(ppiddemo==s)];
            x_dur  = [x_dur ; diseaseDuration(ppiddemo==s)];
            x_sleep= [x_sleep ; psqi_tot(ppidpsqi==s)];
            
            % for final correl w metacog bias
            x_gses = [x_gses  ; gses_tot(ppidgsesmaia==s)];
            x_msse = [x_msse ; msse_tot(ppidgsesmaia==s)];
        else
            disp(['missing some of the data for patient #',num2str(s)])
        end
        % => 66 patients remained in C1 restricted pool
        
    elseif which_hyp == 2
        if ismember(s,ppidhomeo) && ismember(s,ppidmeta2) && ismember(s,ppiddemo) && ismember(s,ppidpsqi)
            % patient s has all required types of data: include in analysis
            
            % define dependent variable
            y_homeo = [y_homeo ; coef_homeo(ppidhomeo==s)];
            
            % define regressors of interest
            x_metacogbias2 = [x_metacogbias2 ; metacogbias2(ppidmeta2==s)];
            x_metacogeff = [x_metacogeff ; metacogeff(ppidmeta2==s)];
            
            % define regressors of no interest
            x_acc2 = [x_acc2 ; acc2(ppidmeta2==s)];
            
            x_age = [x_age ; age(ppiddemo==s)];
            x_sex = [x_sex ; sex(ppiddemo==s)];
            x_medi = [x_medi ; medi(ppiddemo==s)];
            x_medi_immuno = [x_medi_immuno; medi_immuno(ppiddemo==s)];
            x_medi_sedate = [x_medi_sedate; medi_sedate(ppiddemo==s)];
            x_dur  = [x_dur ; diseaseDuration(ppiddemo==s)];
            x_sleep= [x_sleep ; psqi_tot(ppidpsqi==s)];
        else
            disp(['missing some of the data for patient #',num2str(s)])
        end
        % => 61 patients remained in C2 restricted pool
    end
end


if which_hyp == 1
    [rho,pval] = corrcoef(x_gses,x_metacogbias2);
    disp 'Association between GSES and mbias task 2:'
    rho(2)
    pval(2)
    [rho,pval] = corrcoef(x_msse,x_metacogbias2);
    disp 'Association between MS self-efficacy and mbias task 2:'
    rho(2)
    pval(2)
    
    hf = figure('Color','white');
    hold on;
    ylim([10,40]);
    xlim([1,6]);
    plot(x_metacogbias2,x_gses,'ko','LineWidth',1);lsline
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[1,1,1]);
    set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(1,1));
    set(gca,'FontName','Helvetica','FontSize',8);
    ylabel('General self-efficacy score','FontName','Helvetica','FontSize',8);
    xlabel('Metacognitive bias','FontName','Helvetica','FontSize',8);
    set(hf,'PaperPositionMode','manual', ...
        'PaperPosition',[2.5,13,16,4],'PaperUnits','centimeters', ...
        'PaperType','A4','PaperOrientation','portrait');
    figure(hf);
    print('GSES_mbias_correl','-painters','-dpdf');
    
    hf = figure('Color','white');
    hold on;
    ylim([500,1000]);
    xlim([1,6]);
    plot(x_metacogbias2,x_msse,'ko','LineWidth',1);lsline
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[1,1,1]);
    set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(1,1));
    set(gca,'FontName','Helvetica','FontSize',8);
    ylabel('MS self-efficacy score','FontName','Helvetica','FontSize',8);
    xlabel('Metacognitive bias','FontName','Helvetica','FontSize',8);
    set(hf,'PaperPositionMode','manual', ...
        'PaperPosition',[2.5,13,16,4],'PaperUnits','centimeters', ...
        'PaperType','A4','PaperOrientation','portrait');
    figure(hf);
    print('MSSE_mbias_correl','-painters','-dpdf');
end



% -------------- Test our pre-registered hypotheses --------------
% C. Measures of interoception and autonomic
% regulation are related to measures of metacognition

% - Dependent variable: MAIA scores (MAIA subscales 3 [Not-Worrying] and 8
% [Trusting]) (Hypothesis C1) or (Hypothesis C2) first principal component
% of measures of the integrity of homeostatic regulation
% and autonomic function (HRV, deltaBP and deltaHR sitting vs. standing,
% sudomotor activity, sympathetic skin response)

% - Regressors of interest: metacognitive bias (confidence level) and metacognitive efficiency
% For metacognitive bias, we will also add accuracy as a regressor of no interest.

% - Regressors of no interest: duration of disease, age, sex, medication (dose and type).



% zscore our regressors for comparability of regressions coefs:
x_metacogbias2 = (x_metacogbias2-mean(x_metacogbias2))/std(x_metacogbias2);
x_metacogeff = (x_metacogeff-mean(x_metacogeff))/std(x_metacogeff);
x_age = (x_age-mean(x_age))/std(x_age);
x_sex = (x_sex-mean(x_sex))/std(x_sex);
x_medi_immuno = (x_medi_immuno-mean(x_medi_immuno))/std(x_medi_immuno);
x_medi_sedate = (x_medi_sedate-mean(x_medi_sedate))/std(x_medi_sedate);
x_dur  = (x_dur-mean(x_dur))/std(x_dur);
x_sleep  = (x_sleep-mean(x_sleep))/std(x_sleep);

% build regressor matrix
X = [x_metacogbias2 x_metacogeff x_age x_sex x_medi_immuno x_medi_sedate x_dur x_sleep];

if which_hyp == 1
    % perform regression
    [betaC,~,statsC] = glmfit(X,y_maia38);
    fname = './FigureC1';
    % perform GLM to get F-test
    mdl = fitglm(X,y_maia38,'y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8','Distribution','normal');
    p_mdl = coefTest(mdl);
    disp 'Ftest Model C1 with MFIS, pvalue='
    p_mdl
elseif which_hyp == 2
    % perform regression
    [betaC,~,statsC] = glmfit(X,y_homeo);
    fname = './FigureC2';
    % perform GLM to get F-test
    mdl = fitglm(X,y_homeo,'y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8','Distribution','normal');
    p_mdl = coefTest(mdl);
    disp 'Ftest Model C2 with MFIS, pvalue='
    p_mdl
end

disp 'Statistics for each regressor in Model C:'
disp 'betas:'
statsC.beta(2:end)
disp 'tvalues:'
statsC.t(2:end)
disp 'pvalues:'
statsC.p(2:end)



% ---------------- FIGURE C. Metacognition task 2 ----------------

hf = figure('Color','white');

hold on
nbar = length(betaC)-1;
if which_hyp == 1
    ylim([-1,1]);
elseif which_hyp == 2
    ylim([-7,8]);
end
xlim([0.4,nbar+.6]);
pbar = 1/2*(nbar+0.2)/2.2*2;
for i = 1:nbar
    x = betaC(i+1);
    bar(i,mean(x),0.8,'EdgeColor','none','FaceColor',0.5*(rgb(which_hyp,:)+1),'FaceAlpha',0.5);
    pos = i;
    bar(pos,mean(x),0.8,'EdgeColor',rgb(which_hyp,:),'FaceColor','none','LineWidth',1);
    plot(pos*[1,1],[mean(x)-(statsC.se(i+1)) mean(x)+(statsC.se(i+1))],'k','LineWidth',1);
end
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',8);
set(gca,'XTick',1:size(X,2),'XTickLabel',{'m.bias','m.effi','age','sex','mediImmuno','mediSedate','dur','sleep'});
ylabel('Regression coefficients','FontName','Helvetica','FontSize',8);
if which_hyp == 1
    xlabel('Contributors to MAIA scores','FontName','Helvetica','FontSize',8);
    set(gca,'YTick',-1:1,'YTickLabel',{'-1','0','1'});
elseif which_hyp == 2
    xlabel('Contributors to autonomic regulation','FontName','Helvetica','FontSize',8);
end
set(hf,'PaperPositionMode','manual', ...
    'PaperPosition',[2.5,13,16,4],'PaperUnits','centimeters', ...
    'PaperType','A4','PaperOrientation','portrait');
figure(hf);
print(fname,'-painters','-dpdf');

end
