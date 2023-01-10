clear all; close all

% load data file
file = "H:\profiling-resistance-mechanisms\3.resistance-signature\results\performance\otherclones_bortezomibsignature_roc_curve.csv";

% load the table into a variable called "data"
data = readtable(file);

% split data into variables by column
fpr=table2array(data(:,1));
tpr=table2array(data(:,2));
threshold=table2array(data(:,3));
res=table2array(data(:,4));
shuffled=table2array(data(:,5));

% initialize variables
absIX=[];
pointIX=[];

m=1;

% find absolute value of ixazomib thresholds
for i=1:height(data)
    if isequal("ixazomib",res(i)) && isequal("FALSE",shuffled(i))
        absIX(m,1)=i;
        absIX(m,2)=sign(threshold(i))*threshold(i);
        m=m+1;
    end
end

% calculate cut-off point for ixazomib
minVal=min(absIX(:,2));
for j=1:length(absIX)
    if absIX(j,2)==minVal
        pointIX(1,1)=fpr(absIX(j,1));
        pointIX(1,2)=tpr(absIX(j,1));
        IXthresh=threshold(j);
    end
end

% initialize variables
absCB=[];
pointCB=[];

m=1;

% find absolute value of CB5083 thresholds
for i=1:height(data)
    if isequal("cb5083",res(i)) && isequal("FALSE",shuffled(i))
        absCB(m,1)=i;
        absCB(m,2)=sign(threshold(i))*threshold(i);
        m=m+1;
    end
end

% calculate cut-off point for CB5083
minVal=min(absCB(:,2));
for j=1:length(absCB)
    if absCB(j,2)==minVal
        pointCB(1,1)=fpr(absCB(j,1));
        pointCB(1,2)=tpr(absCB(j,1));
        CBthresh=threshold(j);
    end
end

% load data for bortezomib
file = "H:\profiling-resistance-mechanisms\3.resistance-signature\results\performance\otherclones_bortezomibsignature_roc_curve_LAST_BATCH_VALIDATION.csv";

% load the table into a variable called "data"
data = readtable(file);

% split data into columns
fpr=table2array(data(:,1));
tpr=table2array(data(:,2));
threshold=table2array(data(:,3));
shuffled=table2array(data(:,4));

% initialize variables
absBZ=[];
pointBZ=[];

m=1;

% calculate absolute thresholds for bortezomib
for i=1:height(data)
    if isequal("FALSE",shuffled(i))
        absBZ(m,1)=i;
        absBZ(m,2)=sign(threshold(i))*threshold(i);
        m=m+1;
    end
end

% calculate cut-off point for bortezomib
minVal=min(absBZ(:,2));
for j=1:length(absBZ)
    if absBZ(j,2)==minVal
        pointBZ(1,1)=fpr(absBZ(j,1));
        pointBZ(1,2)=tpr(absBZ(j,1));
        BZthresh=threshold(j);
    end
end

% load data for training data set
file = "H:\profiling-resistance-mechanisms\3.resistance-signature\results\performance\bortezomib_roc_curve.csv";

% load the table into a variable called "data"
data = readtable(file);

% separate data into columns
fpr=table2array(data(:,1));
tpr=table2array(data(:,2));
threshold=table2array(data(:,3));
class=table2array(data(:,4));
shuffled=table2array(data(:,5));

% initialize variables
absTRAIN=[];
pointTRAIN=[];
absTEST=[];
pointTEST=[];
absVAL=[];
pointVAL=[];
absHOLD=[];
pointHOLD=[];

m=1;
mm=1;
mmm=1;
mmmm=1;

% calculate absolute thresholds for each dataset
for i=1:height(data)
    if isequal("FALSE",shuffled(i)) && isequal("training",class(i))
        absTRAIN(m,1)=i;
        absTRAIN(m,2)=sign(threshold(i))*threshold(i);
        m=m+1;
    elseif isequal("FALSE",shuffled(i)) && isequal("validation",class(i))
        absVAL(mm,1)=i;
        absVAL(mm,2)=sign(threshold(i))*threshold(i);
        mm=mm+1;
    elseif isequal("FALSE",shuffled(i)) && isequal("test",class(i))
        absTEST(mmm,1)=i;
        absTEST(mmm,2)=sign(threshold(i))*threshold(i);
        mmm=mmm+1;
    elseif isequal("FALSE",shuffled(i)) && isequal("holdout",class(i))
        absHOLD(mmmm,1)=i;
        absHOLD(mmmm,2)=sign(threshold(i))*threshold(i);
        mmmm=mmmm+1;
    end
end

% calculate cut-off points for each dataset
minValTEST=min(absTEST(:,2));
minValTRAIN=min(absTRAIN(:,2));
minValVAL=min(absVAL(:,2));
minValHOLD=min(absHOLD(:,2));

for j=1:length(absTEST)
    if absTEST(j,2)==minValTEST
        pointTEST(1,1)=fpr(absTEST(j,1));
        pointTEST(1,2)=tpr(absTEST(j,1));
        TESTthresh=threshold(j);
    end
end

for j=1:length(absTRAIN)
    if absTRAIN(j,2)==minValTRAIN
        pointTRAIN(1,1)=fpr(absTRAIN(j,1));
        pointTRAIN(1,2)=tpr(absTRAIN(j,1));
        TRAINthresh=threshold(j);
    end
end

for j=1:length(absVAL)
    if absVAL(j,2)==minValVAL
        pointVAL(1,1)=fpr(absVAL(j,1));
        pointVAL(1,2)=tpr(absVAL(j,1));
        VALthresh=threshold(j);
    end
end

for j=1:length(absHOLD)
    if absHOLD(j,2)==minValHOLD
        pointHOLD(1,1)=fpr(absHOLD(j,1));
        pointHOLD(1,2)=tpr(absHOLD(j,1));
        HOLDthresh=threshold(j);
    end
end