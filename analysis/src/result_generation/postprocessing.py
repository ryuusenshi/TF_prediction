import matplotlib.pyplot as plt
import numpy
import pandas

from preprocessing import ORFManipulator


def is_true_tf(true_targets_df, tf_id):
    if (true_targets_df.index==tf_id).any():
        return 1
    else:
        return 0


def is_true_target(true_targets_df, is_tf, tf_id, target_orf_id):
    if not is_tf:
        # not a transcription factor
        return -1
    else:
        try:
            # return 1 if true target, 0 if not
            return true_targets_df.loc[tf_id][str(target_orf_id)] 
        except KeyError:
            # TARGET_ORF_ID not found in ground truth file
            return -2


def create_annotated_analysis_file():
    true_targets_df = pandas.read_csv('../../../data/product/annotated_true_tf_targets.csv', index_col='TF')
    analysis_file = pandas.read_csv('../../../data/product/target_gene_with_header.csv')

    is_tf = analysis_file.apply(lambda x:
        is_true_tf(true_targets_df, x['TF_ORF_ID']),
        axis=1
    )
    analysis_file['IS_TF'] = is_tf

    is_target = analysis_file.apply(lambda x:
        is_true_target(true_targets_df, x['IS_TF'], x['TF_ORF_ID'], x['TARGET_ORF_ID']),
        axis=1
    )
    analysis_file['IS_TARGET'] = is_target

    analysis_file.to_csv('../../../data/product/target_gene_annotated.csv', index=False)
    return analysis_file


def add_header_to_analysis_file():
    target_csv = pandas.read_csv('../../../data/product/target_gene_with_all_p_values.txt', delimiter=' ')
    target_csv.columns = ['TF_ORF_ID', 'SNP_ID', 'SNP_PROBE_NAME', 'TARGET_ORF_ID', 'P_VALUE']
    target_csv.to_csv('../../../data/product/target_gene_with_header.csv', index=False)


# calculate TPR and FPR
def TPR_FPR(P, N, TP, FP):
    return TP/(P+1e-8), FP/(N+1e-8)


def save_plot(name, x, y, title, xlabel, ylabel):
    fig, ax = plt.subplots( nrows=1, ncols=1)
    ax.plot(x, y)
    max_xy = max(x+y)

    # add straight line y=x
    ax.plot([0, max_xy], [0, max_xy], color='k')
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.savefig(name)   # save the figure to file
    plt.close(fig)  
 


###################################
# Marker regression ROCs
#
#

# create ROCs for TFs and SNPs
def create_TARGET_ROCs(df, feature_name, feature_id, tf_id):
    # total positive and negatives for TF
    P, N = len(df[df['IS_TARGET'] == 1]), len(df[df['IS_TARGET'] == 0])

    # characteristic
    p_values = list(df[df['IS_TARGET'] == 1].P_VALUE.values) + [1]
    p_values.sort(reverse=True)

    TPR = []
    FPR = []
    for p_value in p_values:
        df = df[df['P_VALUE'] <= p_value]
        TP = len(df[df['IS_TARGET'] == 1])
        FP = len(df[df['IS_TARGET'] == 0])
 
        TPR_x, FPR_x = TPR_FPR(P, N, TP, FP)
        TPR.append(TPR_x)
        FPR.append(FPR_x)

    basename = '%s_%s_ROC' % (feature_name, feature_id) 
    save_plot('../../../results/%s/%s.png' % (feature_name, basename), FPR, TPR, 'TF target prediction ROC curve for %s_%s' % (feature_name, feature_id), 'FPR', 'TPR')


def create_SNP_TARGET_ROCs():
    df = pandas.read_csv('../../../data/product/target_gene_annotated.csv')

    # only look at TFs
    df = df[df['IS_TF'] == 1]

    for snp_id in df['SNP_ID'].unique():
        snp_df = df[df['SNP_ID'] == snp_id]
        tf_id = snp_df['TF_ORF_ID'].values[0]
 
        # ignore targets that do not exists in ground truth file
        snp_df = snp_df[snp_df['IS_TARGET'] != -2]

        create_TARGET_ROCs(snp_df, 'SNP', snp_id, tf_id)


def create_TF_ROC():
    df = pandas.read_csv('../../../data/product/target_gene_annotated.csv')

    P, N = len(df[df['IS_TF'] == 1]), len(df[df['IS_TF'] == 0])

    p_values = list(numpy.linspace(1, 0.01, num=1000)) + list(numpy.logspace(-2, -18, num=2000))

    TPR = []
    FPR = []
    for p_value in p_values:
        df = df[df['P_VALUE'] <= p_value]
        TP = len(df[df['IS_TF'] == 1])
        FP = len(df[df['IS_TF'] == 0])
 
        TPR_x, FPR_x = TPR_FPR(P, N, TP, FP)
        TPR.append(TPR_x)
        FPR.append(FPR_x)

    save_plot('../../../results/tf_ROC.png', FPR, TPR, 'TF prediction ROC curve', 'FPR', 'TPR')



def create_annotated_snp_to_gene_table():
    orf_manipulator = ORFManipulator()
    true_targets_df = pandas.read_csv('../../../data/product/annotated_true_tf_targets.csv', index_col='TF')
    snp_to_gene_df = pandas.read_csv('../../../data/product/snp_to_gene.csv')

    tf_ids = snp_to_gene_df.apply(lambda x: orf_manipulator.snp_to_orf_id(x['SNP_ID']), axis=1)
    snp_to_gene_df['TF_ORF_ID'] = tf_ids

    is_target = snp_to_gene_df.apply(lambda x: is_true_target(true_targets_df, True, x['TF_ORF_ID'], x['GENE_ID']),
            axis=1
    )
    snp_to_gene_df['IS_TARGET'] = is_target

    # drop targets not appearing in ground truth file
    snp_to_gene_df = snp_to_gene_df[snp_to_gene_df['IS_TARGET'] != -2]
    snp_to_gene_df['total_transitions'] = sum([snp_to_gene_df[x] for x in '0--0;0--1;0--2;1--0;1--1;1--2;2--0;2--1;2--2'.split(';')])
    snp_to_gene_df['major_transitions'] =  sum([snp_to_gene_df[x] for x in '0--0;0--1;0--2'.split(';')])
    snp_to_gene_df['minor_transitions'] =  sum([snp_to_gene_df[x] for x in '2--0;2--1;1--0'.split(';')])
    snp_to_gene_df['major_bias'] = snp_to_gene_df['major_transitions'] / snp_to_gene_df['total_transitions']
    snp_to_gene_df['minor_bias'] = snp_to_gene_df['minor_transitions'] / snp_to_gene_df['total_transitions']

    
    snp_to_gene_df.to_csv('../../../data/product/annotated_snp_to_gene.csv', index=False)
    return snp_to_gene_df


def create_TARGET_ROCs_allele_bias(df, bias_type, feature_name, feature_id, tf_id):
    # total positive and negatives for TF
    P, N = len(df[df['IS_TARGET'] == 1]), len(df[df['IS_TARGET'] == 0])
    
    # characteristic
    biases = numpy.linspace(0, 1, num=1000)


    TPR = []
    FPR = []
    for bias in biases:
        df = df[df[bias_type] >= bias]
        TP = len(df[df['IS_TARGET'] == 1])
        FP = len(df[df['IS_TARGET'] == 0])
 
        TPR_x, FPR_x = TPR_FPR(P, N, TP, FP)
        TPR.append(TPR_x)
        FPR.append(FPR_x)

    basename = '%s_%s_ROC' % (feature_name, feature_id) 
    save_plot('../../../results/BIAS/%s/%s-%s.png' % (feature_name, bias_type, basename), FPR, TPR, 'TF target prediction ROC curve for %s_%s_%s' % (feature_name, bias_type, feature_id), 'FPR', 'TPR')
 

def create_SNP_TARGET_ROCs_allele_bias():
    df = pandas.read_csv('../../../data/product/annotated_snp_to_gene.csv')

    for snp_id in df['SNP_ID'].unique():
        snp_df = df[df['SNP_ID'] == snp_id]
        tf_id = snp_df['TF_ORF_ID'].values[0]
 
        # ignore targets that do not exists in ground truth file
        snp_df = snp_df[snp_df['IS_TARGET'] != -2]

        create_TARGET_ROCs_allele_bias(snp_df, 'major_bias', 'SNP', snp_id, tf_id)
        create_TARGET_ROCs_allele_bias(snp_df, 'minor_bias', 'SNP', snp_id, tf_id)


def sort_importance_table(df):
    # find average importance
    importances = df.apply(lambda x: sum(x)/len(x), axis=0)

    # add same SNPs together
    all_snps = [str(x) for x in importances.index.values if str(x).startswith('SNP')] 
    unique_snps = ['_'.join(x.split('_')[:2]) for x in all_snps]
    for unique_snp in unique_snps:
        copy_snps = [x for x in all_snps if x.startswith(unique_snp)]

        copy_snps = importances.loc[copy_snps]
        importances[unique_snp] = copy_snps.sum() - copy_snps.max()

    importances = importances.sort_values()
    importances.to_csv('../../../results/rf_importances_snps_joined_and_sorted.csv')
    return importances


 
