#- plot_results
import numpy as np
import random
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.metrics import precision_recall_curve, accuracy_score, average_precision_score


def get_dbl_metrics(test, score):
    fig = plt.figure(figsize=(6,3),dpi=100)

    plt.subplot(1,2,1)
    #ROC 
    fpr, tpr, thresholds = metrics.roc_curve(test, score, pos_label=1)
    rauc = metrics.auc(fpr, tpr)
    plt.plot(fpr, tpr,lw=2, label='ROC curve (area = %0.2f)' % rauc)
    plt.plot([0,1],[0,1],'--',color ='black', lw=1)
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.title('ROC (area = %0.2f)' % rauc)

    plt.subplot(1,2,2)
    # precision recall curve
    precision, recall, thresholds = precision_recall_curve(test, score, pos_label=1)
    prauc = metrics.auc(recall, precision)
    plt.plot(recall, precision, lw=2)
    random=len(test[test==1]) / len(test)
    plt.plot([0, 1], [random, random], linestyle='--', label='random', c='black', lw=1)
    plt.xlabel("recall")
    plt.ylabel("precision")
    plt.title("PR curve (area = %0.2f)" % prauc)

    plt.tight_layout()
    
    ap = average_precision_score(test, score)
    
    return rauc, prauc, ap