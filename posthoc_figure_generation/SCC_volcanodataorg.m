clear all; close all

% load signature files
summ="H:\profiling-resistance-mechanisms\3.resistance-signature\results\signatures\signature_summary_bortezomib_signature.csv";
tuk="H:\profiling-resistance-mechanisms\3.resistance-signature\results\signatures\tukey_results_bortezomib_signature.csv";
ccount="H:\profiling-resistance-mechanisms\3.resistance-signature\results\signatures\lm_cell_count_results_bortezomib_signature.csv";
anovadata="H:\profiling-resistance-mechanisms\3.resistance-signature\results\signatures\anova_results_bortezomib_signature.csv";

% load the signature tables into a variables
summary = table2cell(readtable(summ));
tukey=table2cell(readtable(tuk));
cellcount=table2cell(readtable(ccount));
anovatable=readtable(anovadata);

% ID unique features
features=unique(summary(:,2));
terms=unique(tukey(:,2));

% initialize variables excluded from final signature
batch=table;
inclone=table;
status=table;
time=table;
cellcounttable=table;

%initialize variables included in final signature
batchIN=table;
incloneIN=table;
statusIN=table;
timeIN=table;
cellcounttableIN=table;

finalsig=table;

%initialize counting variables
n=1;
nn=1;
nnn=1;
nnnn=1;

nIN=1;
nnIN=1;
nnnIN=1;
nnnnIN=1;

m=1;
mIN=1;
mm=1;

% create table of features included in final signature
for k=1:height(summary)
    if isequal(summary(k,8),{'TRUE'})
        finalsig(mm,1)=summary(k,2);
        mm=mm+1;
    end
end

for i=1:height(tukey)
    % split Batch features into those included vs. excluded in final
    % signature
    if isequal(tukey(i,2),{'Metadata_batch'})
        if ismember(table(tukey(i,8)),finalsig)
            batchIN(nIN,:)=tukey(i,:);
            % replace infinite p values with '20'
            check=table2cell(batchIN(nIN,9));
            if isequal(check,{[Inf]})            
               batchIN(nIN,9)={20};
            end
            nIN=nIN+1;
        else
            batch(n,:)=tukey(i,:);
            % replace infinite p values with '20'
            check=table2cell(batch(n,9));
            if isequal(check,{[Inf]})            
               batch(n,9)={20};
            end
            n=n+1;
        end
    % split Clone ID features into those included vs. excluded in final
    % signature        
    elseif isequal(tukey(i,2),{'Metadata_clone_number'})
        splitclone=split(tukey(i,3),'-');
        % only consider WT-WT comparisons (e.g. WT02 vs WT03), making sure 
        % not to include same clone comparisons (e.g. WT01 vs WT01)
        if contains(splitclone(1), 'WT') && contains(splitclone(2), 'WT') && isequal(splitclone(1),splitclone(2))==0
            if ismember(table(tukey(i,8)),finalsig)
                incloneIN(nnIN,:)=tukey(i,:);
                % replace infinite p values with '20'
                check=table2cell(incloneIN(nnIN,9));
               if isequal(check,{[Inf]})                 
                   incloneIN(nnIN,9)={20};
               end
            nnIN=nnIN+1;
            else
            inclone(nn,:)=tukey(i,:);
            % replace infinite p values with '20'
            check=table2cell(inclone(nn,9));
               if isequal(check,{[Inf]})                 
                   inclone(nn,9)={20};
               end
            nn=nn+1;
            end
        end
    % split Resistance status features into those included vs. excluded in 
    % final signature
    elseif isequal(tukey(i,2),{'Metadata_clone_type_indicator'})
        statRes=tukey(i,3);
        if ismember(table(tukey(i,8)),finalsig)
           statusIN(nnnIN,:)=tukey(i,:);
           % replace infinite p values with '20'
           check=table2cell(statusIN(nnnIN,9));
           if isequal(check,{[Inf]})
             statusIN(nnnIN,9)={20};
           end
        nnnIN=nnnIN+1;
        else
           status(nnn,:)=tukey(i,:);
           % replace infinite p values with '20'
           check=table2cell(status(nnn,9));
           if isequal(check,{[Inf]})
             status(nnn,9)={20};
           end
        nnn=nnn+1;
        end
    % split Treatment Time features into those included vs. excluded in 
    % final signature
    elseif isequal(tukey(i,2),{'Metadata_treatment_time'})
        if ismember(table(tukey(i,8)), finalsig)
           timeIN(nnnnIN,:)=tukey(i,:);
           % replace infinite p values with '20'
           check=table2cell(timeIN(nnnnIN,9));
           if isequal(check,{[Inf]})        
             timeIN(nnnnIN,9)={20};
           end
        nnnnIN=nnnnIN+1;
        else
            time(nnnn,:)=tukey(i,:);
            % replace infinite p values with '20'
            check=table2cell(time(nnnn,9));
            if isequal(check,{[Inf]})        
               time(nnnn,9)={20};
            end
        nnnn=nnnn+1;
        end
    end
end
 
% split Cell Count features into those included vs. excluded in final
% signature
for j=1:height(cellcount)
    if isequal(cellcount(j,2),{'scale(Metadata_cell_count)'})
        if ismember(table(cellcount(j,7)), finalsig)
            cellcounttableIN(mIN,:)=cellcount(j,:);
            mIN=mIN+1;
        else
            cellcounttable(m,:)=cellcount(j,:);
            m=m+1;
        end
    end
end


%% CHECKING VISUALIZATION OF AVERAGE VALUES %%
% 
% incloneFEAT=unique(inclone(:,8));
% incloneINFEAT=unique(incloneIN(:,8));
% 
% batchFEAT=unique(batch(:,8));
% batchINFEAT=unique(batchIN(:,8));
% 
% incloneAVG=table;
% incloneINAVG=table;
% 
% batchAVG=table;
% batchINAVG=table;
% 
% m=1;
% for i=1:height(incloneFEAT)
%     n=1;
%     temp=table;
%     for j=1:height(inclone)
%         if isequal(inclone(j,8), incloneFEAT(i,1))
%             temp(n,1)=inclone(j,4);
%             temp(n,2)=inclone(j,9);
%             n=n+1;
%         end
%     end
%     avgEST=mean(table2array(temp(:,1)));
%     avgP=mean(table2array(temp(:,2)));
%     medianP=median(table2array(temp(:,2)));
%     medianEST=median(table2array(temp(:,1)));
%     incloneAVG(m,1)=incloneFEAT(i,1);
%     incloneAVG(m,2)=array2table(avgEST);
%     incloneAVG(m,3)=array2table(avgP);
%     incloneAVG(m,4)=array2table(medianEST);
%     incloneAVG(m,5)=array2table(medianP);
%     m=m+1;
% end
% 
% m=1;
% for i=1:height(incloneINFEAT)
%     n=1;
%     temp=table;
%     for j=1:height(incloneIN)
%         if isequal(incloneIN(j,8), incloneINFEAT(i,1))
%             temp(n,1)=incloneIN(j,4);
%             temp(n,2)=incloneIN(j,9);
%             n=n+1;
%         end
%     end
%     avgEST=mean(table2array(temp(:,1)));
%     avgP=mean(table2array(temp(:,2)));
%     incloneINAVG(m,1)=incloneINFEAT(i,1);
%     incloneINAVG(m,2)=array2table(avgEST);
%     incloneINAVG(m,3)=array2table(avgP);
%     m=m+1;
% end
% 
% m=1;
% for i=1:height(batchFEAT)
%     n=1;
%     temp=table;
%     for j=1:height(batch)
%         if isequal(batch(j,8), batchFEAT(i,1))
%             temp(n,1)=batch(j,4);
%             temp(n,2)=batch(j,9);
%             n=n+1;
%         end
%     end
%     avgEST=mean(table2array(temp(:,1)));
%     avgP=mean(table2array(temp(:,2)));
%     batchAVG(m,1)=batchFEAT(i,1);
%     batchAVG(m,2)=array2table(avgEST);
%     batchAVG(m,3)=array2table(avgP);
%     m=m+1;
% end
% 
% m=1;
% for i=1:height(batchINFEAT)
%     n=1;
%     temp=table;
%     for j=1:height(batchIN)
%         if isequal(batchIN(j,8), batchINFEAT(i,1))
%             temp(n,1)=batchIN(j,4);
%             temp(n,2)=batchIN(j,9);
%             n=n+1;
%         end
%     end
%     avgEST=mean(table2array(temp(:,1)));
%     avgP=mean(table2array(temp(:,2)));
%     batchINAVG(m,1)=batchINFEAT(i,1);
%     batchINAVG(m,2)=array2table(avgEST);
%     batchINAVG(m,3)=array2table(avgP);
%     m=m+1;
% end
% 
% % anovaFEAT=table;
% % n=1;
% % for i=1:height(anovatable)
% %     if table2array(anovatable(i,9))>4.19423674872383
% %         anovaFEAT(n,1)=anovatable(i,8);
% %         n=n+1;
% %     end
% % end
% % anovaUNIQ=unique(anovaFEAT);