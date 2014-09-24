%% Endogenous attention: Sigmoid Curves

clear all;clc
obs = 'dg';
cont_sub = [.1016 .7367];
nblock_endo =  1;

%% Fixed Parameters
contrast_levels = 2; %c50 and rmax
cue_cond = 2; %valid, invalid
count_block = 0;
%% Data directory
directory = ['/e/4.3/p2/dugue/Documents/MATLAB/TMS/TMS_and_attention/Data/' obs '/TMS1/'];
endo_dir = [directory '/endo_tms_' obs '_b'];

%% Block loop
for block = nblock_endo
    count_block = count_block + 1;
    %% Load data block by block
    if block < 10
        load([endo_dir '0' num2str(block) '.mat'])
    elseif block >= 10
        load([endo_dir num2str(block) '.mat'])
    end
    
    
    %%
    counter = 0;
    for eye = 1:size(real_trial,2)
        if real_trial(1,eye).trialDone == 1
            counter = counter +1;
            eyeOK(counter) = 1;
            eyeOK_idx(eye) = eye;
        end
    end
    eyeOK_idx = find(eyeOK_idx>0);
    %% 
    ntrials = size(trial,2); %determine the number of trials per block
    hit_valid = zeros(1,contrast_levels);
    fa_valid = zeros(1,contrast_levels);
    hit_invalid = zeros(1,contrast_levels);
    fa_invalid = zeros(1,contrast_levels);
    for n = 1:ntrials
        %% Number of trials per cue condition
        a = 0;
        b = 0;
        for nt = 1:ntrials
            if real_trial(1,nt).trialDone == 1
                if (trial(1,nt).cueCond==1 && trial(1,nt).responsecueCond==1) || (trial(1,nt).cueCond==2 && trial(1,nt).responsecueCond==2)...
                        || (trial(1,nt).cueCond==3 && trial(1,nt).responsecueCond==3) || (trial(1,nt).cueCond==4 && trial(1,nt).responsecueCond==4)
                    a = a + 1;
                else
                    b = b + 1;
                end;
            end;
        end;
        
        if n == 1
            %% Separate data into the number of contrast level
            perf_mat_valid = nan(a/contrast_levels,contrast_levels); %create a performance matrix per contrast condition
            perf_mat_invalid = nan(b/contrast_levels,contrast_levels); %create a performance matrix per contrast condition
        end;
        
        %% Separate data into the number of contrast level: HIT vs FALSE ALARME
        % Hit: gabor tilted to the left and subject answers left
        % False alarm: gabor tilted to the right and subject answers left
        hit_a = 0;
        hit_b = 0;
        fa_a = 0;
        fa_b = 0;
        
        for nt = 1:ntrials
            if trial(1,nt).responsecueCond==1 || trial(1,nt).responsecueCond==2
                tmpor = trial(1,nt).oriTest;
            elseif trial(1,nt).responsecueCond==3 || trial(1,nt).responsecueCond==4
                tmpor = trial(1,nt).oriStd;
            end;
            
            if tmpor < 0
                orientation = 1;
            elseif tmpor > 0
                orientation = 2;
            end;
            
            valid = (trial(1,nt).cueCond==1 && trial(1,nt).responsecueCond==1) || (trial(1,nt).cueCond==2 && trial(1,nt).responsecueCond==2)...
                || (trial(1,nt).cueCond==3 && trial(1,nt).responsecueCond==3) || (trial(1,nt).cueCond==4 && trial(1,nt).responsecueCond==4);
            
            invalid = (trial(1,nt).cueCond==1 && trial(1,nt).responsecueCond==3) || (trial(1,nt).cueCond==2 && trial(1,nt).responsecueCond==3)...
                || (trial(1,nt).cueCond==3 && trial(1,nt).responsecueCond==1) || (trial(1,nt).cueCond==4 && trial(1,nt).responsecueCond==1);
            
            if valid
                if orientation == 1
                    hit_a = hit_a + 1;
                end;
                if orientation == 2
                    fa_a = fa_a + 1;
                end
            elseif invalid
                if orientation == 1
                    hit_b = hit_b + 1;
                end
                if orientation == 2
                    fa_b = fa_b + 1;
                end
            end;
        end;
        
        %% Determine the contrast of a given trial
        
        if trial(1,n).stdContrast == cont_sub(1)
            cont = 1;
        elseif trial(1,n).stdContrast == cont_sub(2)
            cont = 2;
        end;
        
        %%
        if trial(1,n).responsecueCond==1 || trial(1,nt).responsecueCond==2
            tmpor = trial(1,n).oriTest;
        elseif trial(1,n).responsecueCond==3 || trial(1,nt).responsecueCond==4
            tmpor = trial(1,n).oriStd;
        end;
        
        if tmpor < 0
            orientation = 1;
        elseif tmpor > 0
            orientation = 2;
        end;
        
        % Hit: gabor tilted to the left and subject answers left
        % False alarm: gabor tilted to the left and subject answers right
        valid = (trial(1,n).cueCond==1 && trial(1,n).responsecueCond==1) || (trial(1,n).cueCond==2 && trial(1,n).responsecueCond==2)...
            || (trial(1,n).cueCond==3 && trial(1,n).responsecueCond==3) || (trial(1,n).cueCond==4 && trial(1,n).responsecueCond==4);
        
        invalid = (trial(1,n).cueCond==1 && trial(1,n).responsecueCond==3) || (trial(1,nt).cueCond==2 && trial(1,n).responsecueCond==3)...
            || (trial(1,n).cueCond==3 && trial(1,n).responsecueCond==1) || (trial(1,nt).cueCond==4 && trial(1,n).responsecueCond==1);
        
        hitValid = 0;
        hitInvalid = 0;
        faValid = 0;
        faInvalid = 0;
        if ~isempty(trial(1,n).responseKey)
            if valid
                if (orientation == 1) && (trial(1,n).responseKey == 1 || trial(1,n).responseKey == 3)
                    hitValid = 1;
                elseif (orientation == 1) && (trial(1,n).responseKey == 2 || trial(1,n).responseKey == 4)
                    faValid = 1;
                end;
            elseif invalid
                if (orientation == 1) && (trial(1,n).responseKey == 1 || trial(1,n).responseKey == 3)
                    hitInvalid = 1;
                elseif (orientation == 1) && (trial(1,n).responseKey == 2 || trial(1,n).responseKey == 4)
                    faInvalid = 1;
                end;
            end
        end
        
        hit_valid(cont) = hit_valid(cont) + hitValid;
        hit_invalid(cont) = hit_invalid(cont) + hitInvalid;
        fa_valid(cont) = fa_valid(cont) + faValid;
        fa_invalid(cont) = fa_invalid(cont) + faInvalid;
        
        %% Determine correct or incorrect for the given trial
        
        if (trial(1,n).cueCond==1 && trial(1,n).responsecueCond==1) || (trial(1,n).cueCond==2 && trial(1,n).responsecueCond==2)...
                        || (trial(1,n).cueCond==3 && trial(1,n).responsecueCond==3) || (trial(1,n).cueCond==4 && trial(1,n).responsecueCond==4)

            if real_trial(1,n).trialDone == 1
                if ~isempty(trial(1,n).cor) || ~isempty(trial(1,n).testChosen)
                    if trial(1,n).testChosen == 1 && trial(1,n).cor == 1
                        perf_valid = 1;
                    else
                        perf_valid = 0;
                    end;
                else
                    perf_valid = 0;
                end;
            else
                perf_valid = 2;
            end
            for a = 1:size(perf_mat_valid,1)
                if isnan(perf_mat_valid(a,cont))
                    perf_mat_valid(a,cont) = perf_valid;
                    break
                else
                    a = a + 1;
                end
            end;
        else
            
            if real_trial(1,n).trialDone == 1
                if ~isempty(trial(1,n).cor) || ~isempty(trial(1,n).testChosen)
                    
                    if trial(1,n).testChosen == 1 && trial(1,n).cor == 1
                        perf_invalid = 1;
                    else
                        perf_invalid = 0;
                    end;
                else
                    perf_invalid = 0;
                end;
            else
                perf_invalid = 2;
            end
            for a = 1:size(perf_mat_invalid,1)
                if isnan(perf_mat_invalid(a,cont))
                    perf_mat_invalid(a,cont) = perf_invalid;
                    break
                else
                    a = a + 1;
                end
            end
        end;
        
    end
    
    if count_block == 1
        valid_perf = perf_mat_valid;
        invalid_perf = perf_mat_invalid;
        
        valid_hit = hit_valid/(hit_a/contrast_levels);
        invalid_hit = hit_invalid/(hit_b/contrast_levels);
        
        valid_fa = fa_valid/(fa_a/contrast_levels);
        invalid_fa = fa_invalid/(fa_b/contrast_levels);
        
    else
        valid_perf = [valid_perf;perf_mat_valid];
        invalid_perf = [invalid_perf;perf_mat_invalid];
        
        valid_fa = [valid_fa;fa_valid/(fa_a/contrast_levels)];
        invalid_fa = [invalid_fa;fa_invalid/(fa_b/contrast_levels)];
        
        valid_hit = [valid_hit;hit_valid/(hit_a/contrast_levels)];
        invalid_hit = [invalid_hit;hit_invalid/(hit_b/contrast_levels)];
        
    end;
end

%% Separate the trials according to the TMS delay

delay = trialSeq(eyeOK_idx,3);
delay_idx = zeros(size(eyeOK,2)/max(unique(delay)),max(unique(delay)));

for a = 1:10
    delay_idx(:,a) = find(delay==a);
end

%% d-prime computation

hit_vec_endo_valid = mean(valid_hit,1);
hit_vec_endo_invalid = mean(invalid_hit,1);

fa_vec_endo_valid = mean(valid_fa,1);
fa_vec_endo_invalid = mean(invalid_fa,1);

for a=1:contrast_levels
    dprime_valid_endo(a)=norminv([hit_vec_endo_valid(a)/100],0,1)-norminv([fa_vec_endo_valid(a)/100],0,1);
    dprime_invalid_endo(a)=norminv([hit_vec_endo_invalid(a)/100],0,1)-norminv([fa_vec_endo_invalid(a)/100],0,1);
end

%% Sigmoidal data fitting
contLevels = cont_sub;
NumPos_valid = zeros(max(unique(delay)),size(contLevels,2));
OutOfNum_valid = zeros(max(unique(delay)),size(contLevels,2));
NumPos_invalid = zeros(max(unique(delay)),size(contLevels,2));
OutOfNum_invalid = zeros(max(unique(delay)),size(contLevels,2));
for a=1:size(contLevels,2)
    for b=1:10
        NumPos_valid(a,b) = size(find(valid_perf(delay_idx(:,b),a)==1),1);
        keyboard
        OutOfNum_valid(a,b) = size(find(valid_perf(delay_idx(:,b),a)==0|valid_perf(delay_idx(:,b),a)==1),1);
        NumPos_invalid(a,b) = size(find(invalid_perf(delay_idx(:,b),a)==1),1);
        OutOfNum_invalid(a,b) = size(find(invalid_perf(delay_idx(:,b),a)==0|invalid_perf(delay_idx(:,b),a)==1),1);
    end
end

%%
params = [.04 .3 .5 .12];
[curveX_valid, fit_valid, thresh75perc_valid, c50_valid, rmax_valid] = pal_fit(contLevels, NumPos_valid, OutOfNum_valid, params);
[curveX_invalid, fit_invalid, thresh75perc_invalid, c50_invalid, rmax_invalid] = pal_fit(contLevels, NumPos_invalid, OutOfNum_invalid, params);

%%%%%%%%%%%%%%
%%% Figures %%
%%%%%%%%%%%%%%
%%% Fitting

xCurve = 0:0.001:8;
% xCurve = 0:8;
xCurveFunct = linspace(0,.4,length(xCurve));

figure;hold on
plot((NumPos_valid./OutOfNum_valid)*100,'o','MarkerFaceColor',[0 .5 .3],'MarkerSize',5,'Color',[0 .5 .3])
plot((NumPos_invalid./OutOfNum_invalid)*100,'o','MarkerFaceColor',[1 0 0],'MarkerSize',5,'Color',[1 0 0])

plot(xCurve,fit_valid*100,'Color',[0 .5 .3],'LineWidth',1.2)
plot(xCurve,fit_invalid*100,'Color',[1 0 0],'LineWidth',1.2)

plot([0 8],[50 50],'--','LineWidth',1,'Color',[0.5 0.5 0.5])

set(gca,'XTick', 1:7)
set(gca,'XTicklabel',[2 4 8 12 16 24 32])
title(['Contrast-response function (' obs ')'],'FontSize',14)
xlabel('Percent contrast' ,'FontSize',10)
ylabel('Percent correct','FontSize',10)
legend(['Valid (ntrials=' num2str(size(valid_perf,1)) ')'],['Invalid (ntrials=' num2str(size(invalid_perf,1)) ')'],'Location','SouthEast')

ylim([30 100])

% namefig=sprintf(['contrastResp_perf_' obs]);
% print ('-djpeg', '-r500',namefig);

%%%%%%%%%%%%%%
%% Figures %%
%%%%%%%%%%%%%%
% d-primes

% initials = [1,1,.5];
% [parameters_valid_dp,otherFitStats] = fitContrastFunction(initials,contLevels,dprime_valid_endo);
% dp_valid_fit = ContrastFunction(parameters_valid_dp, contLevels);
% 
% [parameters_invalid_dp,otherFitStats] = fitContrastFunction(initials,contLevels,dprime_invalid_endo);
% dp_invalid_fit = ContrastFunction(parameters_invalid_dp, contLevels);
% 
% figure; hold on;
% plot(dprime_valid_endo,'o','LineWidth',2,'MarkerFaceColor',[0 .5 .3],'MarkerSize',5,'Color',[0 .5 .3])
% plot(dprime_invalid_endo,'o','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',5,'Color',[1 0 0])
% 
% plot(1:7,dp_valid_fit,'-','LineWidth',2,'Color',[0 .5 .3])
% plot(1:7,dp_invalid_fit,'-','LineWidth',2,'Color',[1 0 0])
% 
% set(gca,'XTicklabel',[0 2 4 8 12 15 24 32 0])
% title('Endogenous attention condition','FontSize',14)
% xlabel('Percent contrast','FontSize',12)
% ylabel('d-prime','FontSize',12)
% legend(['Valid (ntrials=' num2str(size(valid_perf,1)) ')'],['Invalid (ntrials=' num2str(size(invalid_perf,1)) ')'],'Location','SouthEast')
% 
% ylim([-0.2 1])
% 
% namefig=sprintf(['contrastResp_dprime_' obs]);
% print ('-djpeg', '-r500',namefig);
