% Marion Rouault IMEFA study. 2021-2023.

% This is associated with the following paper:
% Rouault, M. et al. Interoceptive and metacognitive facets of
% fatigue in multiple sclerosis. medRxiv 2023.01.23.23284429 (2023).

% This program looks at the behaviour of MS patients in relation to a
% number of psychological and physiological variables and test our
% pre-registered hypothesis B. Perceived fatigue is related to
% measures of metacognition



function [] = ImefaGroup_AnalysisB()

% test hypothesis B1 or B2
which_hyp = 2;
% print pdfs of figures?
print_fig = false;
% analysis of the demographic characteristics?
plot_demo = 0;
% colors
rgb = [[0,0.2,0.8];[0.8,0,0.2];[0.5,0,0.8]; ...
    [0.5,0.8,0];[0.2,0.5,0];[0.2,0.5,0]];




% -------------- retrieve metacognition data for inter-individual analyses:

load groupSES1

ppidmeta1    = groupSES1(:,1);
acc1         = groupSES1(:,2);
metacogbias1 = groupSES1(:,4);

load groupSES2

ppidmeta2    = groupSES2(:,1);
acc2         = groupSES2(:,2);
metacogbias2 = groupSES2(:,4);
metacogeff   = groupSES2(:,11);
dotDiff2     = groupSES2(:,10);



% -------------- retrieve RedCap data for inter-individual analyses:

table = read_RedCap("/Users/marion/Desktop/PDOC/fatigueZurich/fa_analysis/2IMEFA-InitialVar_DATA_LABELS_2021-06-23_1109.csv");

table = table2array(table) ;

% extract data from Metacognition task 2 (ses2)
ses = double(table(:,2));

% 71 participants have demographic data
ppidstring = table(ses==1,1);
agestring  = table(ses==1,3);
sexstring  = table(ses==1,4);
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
diseaseDuration(diseaseDuration<0) = 0; % one person has mistake in dates



% -------------- retrieve RedCap data for fatigue scores

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


% get HDDM parameters for Metacognition task 1
if which_hyp == 1
    
    load('hddm1.mat')
    
    ord_a1          = hddm1(1,:);
    ord_t1          = hddm1(2,:);
    ord_vintercept1 = hddm1(3,:);
    ord_vdelta1     = hddm1(4,:);
    ord_hddm = sort(ppidmeta1);
    
    % order the output of HDDM parameters
    a1          = zeros(length(ppidmeta1),1);
    t1          = zeros(length(ppidmeta1),1);
    vintercept1 = zeros(length(ppidmeta1),1);
    vdelta1     = zeros(length(ppidmeta1),1);
    j = 1;
    for i = 1:length(ppidmeta1)
        a1(j)          = ord_a1(find(ord_hddm==ppidmeta1(i)));
        t1(j)          = ord_t1(find(ord_hddm==ppidmeta1(i)));
        vintercept1(j) = ord_vintercept1(find(ord_hddm==ppidmeta1(i)));
        vdelta1(j)     = ord_vdelta1(find(ord_hddm==ppidmeta1(i)));
        j=j+1;
    end
end


% find restricted pool of participants who have metacog + demo + mfis data

if which_hyp == 1
    x_acc1 = [];
    x_a1 = [];
    x_t1 = [];
    x_vintercept1 = [];
    x_vdelta1 = [];
    x_metacogbias1 = [];
elseif which_hyp == 2
    x_metacogbias2 = [];
    x_metacogeff = [];
    x_acc2 = [];
    x_dots2 = [];
end
y_mfis = [];
x_age = [];
x_sex = [];
x_medi = [];
x_medi_immuno = [];
x_medi_sedate = [];
x_dur = [];
x_sleep = [];


for s = 1:length(ppiddemo) % the largest sample
    if which_hyp == 1
        if ismember(s,ppidmfis) && ismember(s,ppidmeta1) && ismember(s,ppiddemo) && ismember(s,ppidpsqi)
            % patient s has all required types of data: include in analysis
            
            % define dependent variable
            y_mfis = [y_mfis ; mfis_tot(ppidmfis==s)];
            
            % define regressors of interest
            x_metacogbias1 = [x_metacogbias1 ; metacogbias1(ppidmeta1==s)];
            x_a1 = [x_a1 ; a1(ppidmeta1==s)];
            x_t1 = [x_t1 ; t1(ppidmeta1==s)];
            x_vintercept1 = [x_vintercept1 ; vintercept1(ppidmeta1==s)];
            x_vdelta1 = [x_vdelta1 ; vdelta1(ppidmeta1==s)];
            
            % define regressors of no interest
            x_acc1 = [x_acc1 ; acc1(ppidmeta1==s)];
            
            x_age  = [x_age ; age(ppiddemo==s)];
            x_sex  = [x_sex ; sex(ppiddemo==s)];
            x_medi = [x_medi; medi(ppiddemo==s)];
            x_medi_immuno = [x_medi_immuno; medi_immuno(ppiddemo==s)];
            x_medi_sedate = [x_medi_sedate; medi_sedate(ppiddemo==s)];
            x_dur  = [x_dur ; diseaseDuration(ppiddemo==s)];
            x_sleep= [x_sleep;psqi_tot(ppidpsqi==s)];
        else
            disp(['missing some of the data for patient #',num2str(s)])
        end
        % => 55 patients remained in B1 restricted pool
        
    elseif which_hyp == 2
        if ismember(s,ppidmfis) && ismember(s,ppidmeta2) && ismember(s,ppiddemo) && ismember(s,ppidpsqi)
            % patient s has all required types of data: include in analysis
            
            % define dependent variable
            y_mfis = [y_mfis ; mfis_tot(ppidmfis==s)];
            
            % define regressors of interest
            x_metacogbias2 = [x_metacogbias2 ; metacogbias2(ppidmeta2==s)];
            x_metacogeff = [x_metacogeff ; metacogeff(ppidmeta2==s)];
            
            % define regressors of no interest
            x_acc2 = [x_acc2 ; acc2(ppidmeta2==s)];
            x_dots2 = [x_dots2; dotDiff2(ppidmeta2==s)];
            
            x_age = [x_age ; age(ppiddemo==s)];
            x_sex = [x_sex ; sex(ppiddemo==s)];
            x_medi = [x_medi ; medi(ppiddemo==s)];
            x_medi_immuno = [x_medi_immuno ; medi_immuno(ppiddemo==s)];
            x_medi_sedate = [x_medi_sedate ; medi_sedate(ppiddemo==s)];
            x_dur  = [x_dur ; diseaseDuration(ppiddemo==s)];
            x_sleep= [x_sleep;psqi_tot(ppidpsqi==s)];
        else
            disp(['missing some of the data for patient #',num2str(s)])
        end
        % => 56 patients remained in B2 restricted pool
    end
end







% -------------- Test our pre-registered hypotheses --------------
% B. Perceived fatigue is related to measures of metacognition


if which_hyp == 1
    
    % -------------- TEST B. Metacognition task 1 --------------
    
    % - metacognitive bias is negatively associated with fatigue
    
    % zscore our regressors for comparability of regressions coefs:
    x_metacogbias1 = (x_metacogbias1-mean(x_metacogbias1))/std(x_metacogbias1);
    x_acc1 = (x_acc1-mean(x_acc1))/std(x_acc1);
    x_a1 = (x_a1-mean(x_a1))/std(x_a1);
    x_t1 = (x_t1-mean(x_t1))/std(x_t1);
    x_vintercept1 = (x_vintercept1-mean(x_vintercept1))/std(x_vintercept1);
    x_vdelta1 = (x_vdelta1-mean(x_vdelta1))/std(x_vdelta1);
    x_age = (x_age-mean(x_age))/std(x_age);
    x_sex = (x_sex-mean(x_sex))/std(x_sex);
    x_medi_immuno = (x_medi_immuno-mean(x_medi_immuno))/std(x_medi_immuno);
    x_medi_sedate = (x_medi_sedate-mean(x_medi_sedate))/std(x_medi_sedate);
    x_dur = (x_dur-mean(x_dur))/std(x_dur);
    x_sleep = (x_sleep-mean(x_sleep))/std(x_sleep);
    
    % build regressor matrix
    X1 = [x_metacogbias1 x_acc1 x_age x_sex x_medi_immuno x_medi_sedate x_dur x_sleep];
    X1hddm = [x_metacogbias1 x_a1 x_t1 x_vintercept1 x_vdelta1 x_age x_sex x_medi_immuno x_medi_sedate x_dur x_sleep];
    
    % perform regression
    [betaB1,~,statsB1] = glmfit(X1,y_mfis);
    [betaB1hddm,~,statsB1hddm] = glmfit(X1hddm,y_mfis);
    
    % perform GLM to get F-test
    mdl = fitglm(X1,y_mfis,'y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8','Distribution','normal');
    p_mdl = coefTest(mdl);
    disp 'Ftest Model B1 with MFIS, pvalue='
    p_mdl
    
    mdlhddm = fitglm(X1hddm,y_mfis,'y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11','Distribution','normal');
    p_mdlhddm = coefTest(mdlhddm);
    disp 'Ftest Model B1 HDDM with MFIS, pvalue='
    p_mdlhddm
    
    disp 'Statistics for each regressor in Model B1:'
    disp 'betas:'
    statsB1.beta(2:end)
    disp 'tvalues:'
    statsB1.t(2:end)
    disp 'pvalues:'
    statsB1.p(2:end)
    
    % examine Benjamini-Hochler FDR correction:
    orig_pvalues = [[1:length(statsB1.p(2:end))]' statsB1.p(2:end)];
    % order the pvalues:
    ordered_pvalues = sortrows(orig_pvalues,2);
    % add rank:
    ordered_pvalues = [ordered_pvalues [1:length(statsB1.p(2:end))]'];
    % what percentage of FDR rate you decide
    fdr_rate = .05;
    % how many tests do you want to correct for
    m = 2; % here, accu and metacogbias
    % calculate each individual p-value?s Benjamini-Hochberg critical value
    ordered_pvalues(:,4) = ordered_pvalues(:,3)./m*fdr_rate;
    % identify the highest p-value that is also smaller than the critical value
    ordered_pvalues
    
    
    
    disp 'Statistics for each regressor in Model B1hddm:'
    disp 'betas:'
    statsB1hddm.beta(2:end)
    disp 'tvalues:'
    statsB1hddm.t(2:end)
    disp 'pvalues:'
    statsB1hddm.p(2:end)
    
    % examine Benjamini-Hochler FDR correction:
    orig_pvalues = [[1:length(statsB1hddm.p(2:end))]' statsB1hddm.p(2:end)];
    % order the pvalues:
    ordered_pvalues = sortrows(orig_pvalues,2);
    % add rank:
    ordered_pvalues = [ordered_pvalues [1:length(statsB1hddm.p(2:end))]'];
    % what percentage of FDR rate you decide
    fdr_rate = .05;
    % how many tests do you want to correct for
    m = 5; % here, 4 DDM params and metacogbias
    % calculate each individual p-values Benjamini-Hochberg critical value
    ordered_pvalues(:,4) = ordered_pvalues(:,3)./m*fdr_rate;
    % identify the highest p-value that is also smaller than the critical value
    ordered_pvalues
    
    
    
    % -------------- FIGURE B. Metacognition task 1 --------------
    
    hf = figure('Color','white');
    
    hold on
    nbar = length(betaB1)-1;
    xlim([0.4,nbar+.6]);
    ylim([-11,12]);
    pbar = 1/2*(nbar+0.2)/2.2*2;
    for i = 1:nbar
        x = betaB1(i+1);
        bar(i,mean(x),0.8,'EdgeColor','none','FaceColor',0.5*(rgb(which_hyp,:)+1),'FaceAlpha',0.5);
        pos = i;
        bar(pos,mean(x),0.8,'EdgeColor',rgb(which_hyp,:),'FaceColor','none','LineWidth',1);
        plot(pos*[1,1],[mean(x)-(statsB1.se(i+1)) mean(x)+(statsB1.se(i+1))],'k','LineWidth',1);
    end
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
    set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
    set(gca,'FontName','Helvetica','FontSize',10);
    set(gca,'XTick',1:size(X1,2),'XTickLabel',{'m.bias','accu','age','sex','mediImmuno','mediSedate','duration','sleep'});
    ylabel('Regression coefficients','FontName','Helvetica','FontSize',10);
    xlabel('Contribution to fatigue scores (MFIS)','FontName','Helvetica','FontSize',10);
    
    set(hf,'PaperPositionMode','manual', ...
        'PaperPosition',[2.5,13,16,4],'PaperUnits','centimeters', ...
        'PaperType','A4','PaperOrientation','portrait');
    figure(hf);
    fname = './FigureB1';
    if print_fig
        print(fname,'-painters','-dpdf');
    end
    
    hf = figure('Color','white');
    
    hold on
    nbar = length(betaB1hddm)-1;
    xlim([0.4,nbar+.6]);
    ylim([-11,12]);
    pbar = 1/2*(nbar+0.2)/2.2*2;
    for i = 1:nbar
        x = betaB1hddm(i+1);
        bar(i,mean(x),0.8,'EdgeColor','none','FaceColor',0.5*(rgb(which_hyp,:)+1),'FaceAlpha',0.5);
        pos = i;
        bar(pos,mean(x),0.8,'EdgeColor',rgb(which_hyp,:),'FaceColor','none','LineWidth',1);
        plot(pos*[1,1],[mean(x)-(statsB1hddm.se(i+1)) mean(x)+(statsB1hddm.se(i+1))],'k','LineWidth',1);
    end
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
    set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
    set(gca,'FontName','Helvetica','FontSize',10);
    set(gca,'XTick',1:size(X1hddm,2),'XTickLabel',{'m.bias','a','t','vintercept','vdelta','age','sex','mediImmuno','mediSedate','duration','sleep'});
    ylabel('Regression coefficients','FontName','Helvetica','FontSize',10);
    xlabel('Contribution to fatigue scores (MFIS)','FontName','Helvetica','FontSize',10);
    
    set(hf,'PaperPositionMode','manual', ...
        'PaperPosition',[2.5,13,16,4],'PaperUnits','centimeters', ...
        'PaperType','A4','PaperOrientation','portrait');
    figure(hf);
    fname = './FigureB1hddm';
    if print_fig
        print(fname,'-painters','-dpdf');
    end
    
    
elseif which_hyp == 2
    
    % -------------- TEST B. Metacognition task 2 --------------
    
    % - Dependent variable: MFIS scores
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
    x_dur = (x_dur-mean(x_dur))/std(x_dur);
    x_sleep = (x_sleep-mean(x_sleep))/std(x_sleep);
    
    % build regressor matrix
    X2 = [x_metacogbias2 x_metacogeff x_age x_sex x_medi_immuno x_medi_sedate x_dur x_sleep];
    
    % perform regression
    [betaB2,~,statsB2] = glmfit(X2,y_mfis);
    
    % perform GLM to get F-test
    mdl = fitglm(X2,y_mfis,'y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8','Distribution','normal');
    p_mdl = coefTest(mdl);
    disp 'Ftest Model B2 with MFIS, pvalue='
    p_mdl
    
    disp 'Statistics for each regressor in Model B2:'
    disp 'betas:'
    statsB2.beta(2:end)
    disp 'tvalues:'
    statsB2.t(2:end)
    disp 'pvalues:'
    statsB2.p(2:end)
    
    
    % -------------- FIGURE B. Metacognition task 2 --------------
    
    hf = figure('Color','white');
    
    hold on
    nbar = length(betaB2)-1;
    ylim([-9,12]);
    xlim([0.4,nbar+.6]);
    pbar = 1/2*(nbar+0.2)/2.2*2;
    for i = 1:nbar
        x = betaB2(i+1);
        bar(i,mean(x),0.8,'EdgeColor','none','FaceColor',0.5*(rgb(which_hyp,:)+1),'FaceAlpha',0.5);
        pos = i;
        bar(pos,mean(x),0.8,'EdgeColor',rgb(which_hyp,:),'FaceColor','none','LineWidth',1);
        plot(pos*[1,1],[mean(x)-(statsB2.se(i+1)) mean(x)+(statsB2.se(i+1))],'k','LineWidth',1);
    end
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
    set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
    set(gca,'FontName','Helvetica','FontSize',10);
    set(gca,'XTick',1:size(X2,2),'XTickLabel',{'m.bias','m.effi','age','sex','mediImmuno','mediSedate','duration','sleep'});
    ylabel('Regression coefficients','FontName','Helvetica','FontSize',10);
    xlabel('Contributors to fatigue scores (MFIS)','FontName','Helvetica','FontSize',10);
    
    
    set(hf,'PaperPositionMode','manual', ...
        'PaperPosition',[2.5,13,16,4],'PaperUnits','centimeters', ...
        'PaperType','A4','PaperOrientation','portrait');
    figure(hf);
    fname = './FigureB2';
    if print_fig
        print(fname,'-painters','-dpdf');
    end
    
end



% ---------------- plot a number of sanity checks ----------------

if plot_demo
    if which_hyp == 1
        figure;
        plot(x_age,x_acc1,'bo');lsline
        ylabel('Performance Metacog task1','FontName','Helvetica','FontSize',15);
        xlabel('Patient age','FontName','Helvetica','FontSize',15);
        set(gca,'FontName','Helvetica','FontSize',15);
        
    elseif which_hyp == 2
        figure;
        
        subplot(2,2,1)
        histfit(y_mfis,12)
        xlabel('MFIS total score','FontName','Helvetica','FontSize',15);
        ylabel('Patient frequency','FontName','Helvetica','FontSize',15);
        set(gca,'FontName','Helvetica','FontSize',15);
        
        subplot(2,2,2)
        plot(x_age,y_mfis,'bo');lsline
        ylabel('MFIS total score','FontName','Helvetica','FontSize',15);
        xlabel('Patient age','FontName','Helvetica','FontSize',15);
        set(gca,'FontName','Helvetica','FontSize',15);
        
        ttt = .17;
        subplot(2,2,3)
        hold on;
        errorbar(1+ttt:2+ttt,[mean(x_metacogbias2(x_sex == 2)) mean(x_metacogbias2(x_sex == 1))], ...
            [std(x_metacogbias2(x_sex == 2))/sqrt(length(x_metacogbias2(x_sex == 2))) ...
            std(x_metacogbias2(x_sex == 1))/sqrt(length(x_metacogbias2(x_sex == 1)))], ...
            'LineWidth',2)
        plot(3-x_sex,x_metacogbias2,'co')
        ylabel('Metacognitive bias','FontName','Helvetica','FontSize',15);
        set(gca,'FontName','Helvetica','FontSize',15);
        set(gca,'XTick',1:2,'XTickLabel',{'Female','Male'},'FontSize',15);
        axis([.6 2.4 1 6])
        hold off
    end
    
    
    % Before analyses, characterisation of demographics of the sample:
    
    hf = figure('Color','white');
    hold on
    histfit(age);
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[1,1,1]);
    set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(1,1));
    set(gca,'FontName','Helvetica','FontSize',8);
    ylabel('Number of persons','FontSize',8);
    xlabel('Age','FontSize',8);
    set(hf,'PaperPositionMode','manual', ...
        'PaperPosition',[2.5,13,16,4],'PaperUnits','centimeters', ...
        'PaperType','A4','PaperOrientation','portrait');
    figure(hf);
    fname = './FigureDemoAge';
    print(fname,'-painters','-dpdf');
    
    
    hf = figure('Color','white');
    hold on
    nbar = 2;
    ylim([0,62]);
    xlim([0.4,nbar+.6]);
    pbar = 1/2*(nbar+0.2)/2.2*2;
    for i = 1:nbar
        x = length(find(sex==i));
        bar(i,mean(x),0.8,'EdgeColor','none','FaceColor',0.5*(rgb(2+i,:)+1),'FaceAlpha',0.5);
        bar(i,mean(x),0.8,'EdgeColor',rgb(2+i,:),'FaceColor','none','LineWidth',1);
    end
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
    set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
    set(gca,'FontName','Helvetica','FontSize',8);
    ylabel('Number of persons','FontSize',8);
    xlabel(' ','FontSize',8);
    set(gca,'XTick',1:2,'XTickLabel',{'men','women'});
    set(hf,'PaperPositionMode','manual', ...
        'PaperPosition',[2.5,13,16,4],'PaperUnits','centimeters', ...
        'PaperType','A4','PaperOrientation','portrait');
    figure(hf);
    fname = './FigureDemoSex';
    print(fname,'-painters','-dpdf');
    
    
    hf = figure('Color','white');
    hold on
    nbar = 2;
    ylim([0,62]);
    xlim([0.4,nbar+.6]);
    pbar = 1/2*(nbar+0.2)/2.2*2;
    for i = 1:nbar
        x = length(find(medi_immuno==(i-1)));
        bar(i,mean(x),0.8,'EdgeColor','none','FaceColor',0.5*(rgb(5,:)+1),'FaceAlpha',0.5);
        bar(i,mean(x),0.8,'EdgeColor',rgb(5,:),'FaceColor','none','LineWidth',1);
    end
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
    set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
    set(gca,'FontName','Helvetica','FontSize',8);
    ylabel('Number of persons','FontSize',8);
    xlabel('   ','FontSize',8);
    set(gca,'XTick',1:2,'XTickLabel',{'not medic.','medicated'});
    set(hf,'PaperPositionMode','manual', ...
        'PaperPosition',[2.5,13,16,4],'PaperUnits','centimeters', ...
        'PaperType','A4','PaperOrientation','portrait');
    figure(hf);
    fname = './FigureDemoMedicationImmuno';
    print(fname,'-painters','-dpdf');
    
    
    hf = figure('Color','white');
    hold on
    nbar = 2;
    ylim([0,62]);
    xlim([0.4,nbar+.6]);
    pbar = 1/2*(nbar+0.2)/2.2*2;
    for i = 1:nbar
        x = length(find(medi_sedate==(i-1)));
        bar(i,mean(x),0.8,'EdgeColor','none','FaceColor',0.5*(rgb(5,:)+1),'FaceAlpha',0.5);
        bar(i,mean(x),0.8,'EdgeColor',rgb(5,:),'FaceColor','none','LineWidth',1);
    end
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
    set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
    set(gca,'FontName','Helvetica','FontSize',8);
    ylabel('Number of persons','FontSize',8);
    xlabel('   ','FontSize',8);
    set(gca,'XTick',1:2,'XTickLabel',{'not medic.','medicated'});
    set(hf,'PaperPositionMode','manual', ...
        'PaperPosition',[2.5,13,16,4],'PaperUnits','centimeters', ...
        'PaperType','A4','PaperOrientation','portrait');
    figure(hf);
    fname = './FigureDemoMedicationSedate';
    print(fname,'-painters','-dpdf');
    
    
    hf = figure('Color','white');
    hold on
    xlim([0,35]);
    hist(years(consent_day-diagnostic_day));
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[1,1,1]);
    set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(1,1));
    set(gca,'FontName','Helvetica','FontSize',8);
    ylabel('Number of persons','FontSize',8);
    xlabel('Disease duration (years)','FontSize',8);
    set(hf,'PaperPositionMode','manual', ...
        'PaperPosition',[2.5,13,16,4],'PaperUnits','centimeters', ...
        'PaperType','A4','PaperOrientation','portrait');
    figure(hf);
    fname = './FigureDiseaseDuration';
    print(fname,'-painters','-dpdf');
end

end
