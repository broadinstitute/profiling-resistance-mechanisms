AUROC_README.md

AUROC values calculated using roc_auc_score were different from values calculated using 
Graphpad PRISM when double checking data for figure generation. This is a known issue 
for roc_auc_score (https://github.com/scikit-learn/scikit-learn/issues/10914). The values 
from PRISM more closely corresponded with the expected values based on qualitative 
analysis of the ROC curves and were therefore used for publication.