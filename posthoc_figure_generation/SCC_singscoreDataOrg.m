clear all; close all

% load file of interest (uncomment as necessary)
%file = "H:\profiling-resistance-mechanisms\3.resistance-signature\results\singscore\singscore_results_otherclones_UNIQUE.csv";
%file = "H:\profiling-resistance-mechanisms\3.resistance-signature\results\singscore\singscore_results_LAST_BATCH_VALIDATIONotherclones_UNIQUE.csv";
file="H:\profiling-resistance-mechanisms\3.resistance-signature\results\singscore\singscore_resultsbortezomib";

% load the table into a variable called "data"
data = readtable(file);


% extract columns from the data table
col2=table2array(data(:,2)); % plate number
col20=table2array(data(:,20)); % singscore
col8=table2array(data(:,8)); % clone ID
col18=table2array(data(:,18)); % model

% create an array containing the unique entries
%uniq2=unique(col2);
uniq8 = unique(col8);
uniq18=unique(col18);

%initialize variables

% CBtable=table;
% IXtable=table;
% BZtable=table;

holdout=table;
test=table;
train=table;
validate=table;
inference=table;

%% CODE FOR OTHER CLONES AND LAST BATCH VALIDATION DATA %%
% % uncomment as necessary
%
% CB=[218696, 218774, 218852, 218856];
% IX=[218698, 218854, 218858];
% 
% % split singscores by clone ID for CB5083
% for k=1:height(uniq8)
%     n=2;
%     for j=1:height(data)
%         if isequal(col8(j),uniq8(k)) && ismember(col2(j), CB)
%             CBtable{k,1}=strrep(col8(j), "'",'');
%             CBtable{k,n}=col20(j);
%             n=n+1;
%         end
%     end
% end
%
% % split singscores by clone ID for ixazomib
% for k=1:height(uniq8)
%     n=2;
%     for j=1:height(data)
%         if isequal(col8(j),uniq8(k)) && ismember(col2(j), IX)
%             IXtable{k,1}=strrep(col8(j),"'",'');
%             IXtable{k,n}=col20(j);
%             n=n+1;
%         end
%     end
% end
% 
% % split singscores by clone ID for bortezomib
% for k=1:height(uniq8)
%     n=2;
%     for j=1:height(data)
%         if isequal(col8(j),uniq8(k))
%             BZtable{k,1}=strrep(col8(j),"'",'');
%             BZtable{k,n}=col20(j);
%             n=n+1;
%         end
%     end
% end
% 

% comment as necessary
% split singscores by model and clone ID for training dataset
for k=1:height(uniq8)
    n=2;
    nn=2;
    nnn=2;
    nnnn=2;
    for j=1:height(data)
        if isequal(col8(j),uniq8(k)) && isequal({'holdout'}, col18(j))
            holdout{k,1}=strrep(col8(j),"'",'');
            holdout{k,n}=col20(j);
            n=n+1;
        elseif isequal(col8(j),uniq8(k)) && isequal({'test'}, col18(j))
            test{k,1}=strrep(col8(j),"'",'');
            test{k,nn}=col20(j);
            nn=nn+1;
        elseif isequal(col8(j),uniq8(k)) && isequal({'training'}, col18(j))
            train{k,1}=strrep(col8(j),"'",'');
            train{k,nnn}=col20(j);
            nnn=nnn+1;
        elseif isequal(col8(j),uniq8(k)) && isequal({'validation'}, col18(j))
            validate{k,1}=strrep(col8(j),"'",'');
            validate{k,nnnn}=col20(j);
            nnnn=nnnn+1; 
        end
    end
end


