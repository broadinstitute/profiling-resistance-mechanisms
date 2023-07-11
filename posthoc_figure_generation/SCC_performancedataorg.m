clear all; close all
% tsv files converted to csv for analysis

% load data file
file="H:\profiling-resistance-mechanisms\3.resistance-signature\results\performance\total_bortezomib_shuffle_metric_performance.csv";

% load the table into a variable called "data"
data = readtable(file);

% split data into arrays
model=table2array(data(:,1));
ap=table2array(data(:,2));
metric=data(:,3);

% generate lists of unique entries in model and ap
modelset=unique(model);
apset=unique(ap);

%initialize variables
accT=table;
preT=table;

%split data into tables of accuracy and precision for each unique dataset
for j=1:length(modelset)
    n=2;
    m=2;
    for k=1:height(metric)
        if isequal(ap(k),{'accuracy'}) && isequal(model(k),modelset(j))
            accT(j,1)=model(k,1);
            accT(j,n)=metric(k,1);
            n=n+1;
        elseif isequal(ap(k),{'avg_precision'}) && isequal(model(k),modelset(j))
            preT(j,1)=model(k,1);
            preT(j,m)=metric(k,1);
            m=m+1;
        end
    end
end