import csv
import pandas

from collections import defaultdict


class TFManipulator(object):
    def __init__(self):
        orf_manipulator = ORFManipulator()
        tfs_with_snps = pandas.read_csv('../../../data/product/tf_snp_count.csv')
        tfs_with_snps = tfs_with_snps[tfs_with_snps['SNP_COUNT']>0]

        snps = []
        for tf_id in tfs_with_snps['ID']:
            snps = snps + list(orf_manipulator.orf_id_to_snp(tf_id).values)

        self.tfs_with_snps = tfs_with_snps
        self.snps = snps


class TargetManipulator(object):
    def __init__(self):
        self.tf_targets = pandas.read_csv('../../../data/product/annotated_true_tf_targets.csv', index_col='TF')


class ORFManipulator(object):
    _orf_to_snp_df = None

    def __init__(self):
        if self._orf_to_snp_df is None:
            try:
                self._orf_to_snp_df = pandas.read_csv('../../../data/product/snp_orf.csv')
            except:
                raise Warning("Make sure to run create_snp_orf_file before using ORFManipulator")
            self._orf_info = pandas.read_csv('../../../data/id_to_gene_name.csv')

    def gene_name_to_orf_id(self, gene_name):
        if (self._orf_info['Platform_ORF']==gene_name.upper()).any():
            src = 'Platform_ORF'
        elif (self._orf_info['Gene title']==gene_name).any():
            src = 'Gene title'
        else:
            src = 'Gene symbol'

        df = self._orf_info.loc[self._orf_info[src] == gene_name.upper()]
 
        try:
            return df['ID'].values[0]
        except IndexError:
            return -1

    def orf_id_to_snp(self, orf_id):
        df = self._orf_to_snp_df.loc[self._orf_to_snp_df['ID'] == int(orf_id)]

        return df['SNP_ID']

    def snp_to_orf_id(self, snp_id):
        df = self._orf_to_snp_df.loc[self._orf_to_snp_df['SNP_ID'] == int(snp_id)]

        return df['ID'].values[0]

    def orf_id_to_orf_name(self, orf_id):
        df = self._orf_info.loc[self._orf_info['ID'] == int(orf_id)]

        return df['Platform_ORF'].values[0]

    def orf_name_to_orf_id(self, orf_name):
        df = self._orf_info.loc[self._orf_info['Platform_ORF'] == orf_name]

        return df['ID'].values[0]

    def orf_id_to_gene_name(self, orf_id):
        df = self._orf_info.loc[self._orf_info['ID'] == int(orf_id)]

        try:
            return df['Gene symbol'].values[0]
        except:
            return None



def create_snp_orf_file():
    # inputs
    orf_df = pandas.read_csv('../../../data/expression_location', delimiter='\t')
    snp_df = pandas.read_csv('../../../data/snp_position.txt', delimiter='\t')

    # output
    f = open('../../../data/product/snp_orf.csv', 'w')
    writer = csv.writer(f, delimiter=',')
    writer.writerow(['SNP_ID', 'ID'])

    for snp_id, snp in snp_df.iterrows():
        orf = orf_df.loc[
                orf_df['Chr'] == snp['chromosome']
            ].loc[
                snp['position'] >= orf_df['Start']
            ].loc[
                snp['position'] <= orf_df['Stop']
            ]
        if len(orf['ID']) == 0:
            orf = orf_df.loc[
                orf_df['Chr'] == snp['chromosome']
            ].loc[
                snp['position'] <= orf_df['Start']
        ].loc[
            snp['position'] >= orf_df['Stop']
        ]
        
        if len(orf['ID']) > 0:
            for orf_id in orf['ID']:
                writer.writerow([snp_id, orf_id])
    f.close()
    return pandas.read_csv('../../../data/product/snp_orf.csv')


def create_annotated_true_tf_target_table():
    orf_manipulator = ORFManipulator()
    # input
    df = pandas.read_csv('../../../data/RegulationMatrix.csv', delimiter=";")
    row_id_column = df.columns.values[0]

    TFs = df[row_id_column].apply(lambda x: x[:-1] if x.endswith('p') else x)
    TFids = TFs.apply(lambda x: int(orf_manipulator.gene_name_to_orf_id(x)))

    df[row_id_column] = TFids
    df = df.rename(index=int, columns={row_id_column: 'TF'})
    df = df.set_index('TF')

    df.columns = df.columns.map(lambda x: orf_manipulator.gene_name_to_orf_id(x))
    # drop rows - unidentified TFs
    df = df.drop(-1, axis=0)
    # drop columns - unidentified target genes
    df = df.drop(-1, axis=1)


    # output
    df.to_csv('../../../data/product/annotated_true_tf_targets.csv')
    return df




# ALLELE TRANSITION BIAS
log_threshold = 2


def create_upregulated_targets_file():
    # input file
    exp_df = pandas.read_csv('../../../data/expression_matrix.txt', delimiter='\t')

    # output file
    out_filename = '../../../data/product/upregulated_targets.csv'
    f = open(out_filename, 'w')
    w = csv.writer(f)
    # write header to output file
    w.writerow(['ID', 'DOWN_CONDITION--UP_CONDITION'])

    conditions = exp_df.columns.values[1:]
    n_cond = len(conditions)
    for index, row in exp_df.iterrows():
        upregulated_transitions = []
        for i in range(n_cond):
            for j in range(i+1, n_cond):
                cond_i = conditions[i]
                cond_j = conditions[j]
                exp_i = row[cond_i]
                exp_j = row[cond_j]

                if exp_i-exp_j >= log_threshold:
                    cond_transition = '%s--%s' % (cond_j, cond_i)
                elif exp_j-exp_i >= log_threshold:
                    cond_transition = '%s--%s' % (cond_i, cond_j)
                else:
                    cond_transition = None

                if cond_transition is not None:
                    upregulated_transitions.append(cond_transition)
        w.writerow([
            '%d' % row['ID'],
            ';'.join(upregulated_transitions)
        ])
    f.close()
    return pandas.read_csv('../../../data/product/upregulated_targets.csv')


def create_snp_to_gene_file():
    # -- is shorthand for 00, 11, and 22 (no change)
    # ^^ is shorhand for 02, 01, 12 (minor to major)
    # vv is shorthand for 20, 10, 21 (major to minor)
    # IMPORTANT: SNP IDS are counted from 0
    orf_manipulator = ORFManipulator()
    def get_up_transitions(upregulated_series):
        transition_list = upregulated_series.loc['DOWN_CONDITION--UP_CONDITION']
        return set(transition_list.split(';'))

    def multipleReplace(text, wordDict):
        for key in wordDict:
            text = text.replace(key, str(wordDict[key]))
        return text

    def allele_transition_counts(text):
        allele_transitions = ['0--2', '0--1', '1--2', '2--0', '1--0', '2--1', '0--0', '1--1', '2--2']

        return dict([(x, text.count(x)) for x in allele_transitions])

    filename = '../../../data/product/snp_to_gene.csv'
    f = open(filename, 'w')
    fieldnames = ['SNP_ID', 'GENE_ID', '0--2', '0--1', '1--2', '2--0', '1--0', '2--1', '0--0', '1--1', '2--2'] 
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()


    # only SNPs for TFs
    snps = snp_ids_for_tfs()
    snp_df = pandas.read_csv('../../../data/matrix_genotypes.txt', delimiter='\t')

    snp_df = snp_df.loc[snps]
    conditions = snp_df.columns.values[1:]
    # add SNP_ID to the snp_df dataframe
    snps = pandas.Series(snps, name='SNP_ID', index=snp_df.index)
    snp_df = snp_df.join(snps)


    upregulated_df = pandas.read_csv('../../../data/product/upregulated_targets.csv')
    upregulated_df = upregulated_df.dropna()
    for _, snp in snp_df.iterrows(): 
        snp_name = snp['SNP_ID']
        snp_values = snp.to_dict()
        print("Working on SNP: %d" % snp_name)
        snp_gene_id = orf_manipulator.snp_to_orf_id(snp_name)

        snp_dict = {}
        gene_ids = set(upregulated_df['ID']) 
        for gene_id in gene_ids:
            snp_dict[gene_id] = defaultdict(int)
            snp_dict[gene_id].update({
                'SNP_ID': snp_name,
                'GENE_ID': gene_id,
            })
     
        # skip snps whose genes showed upregulation
        try:
            snp_gene_upregulated_transitions = get_up_transitions(upregulated_df.loc[snp_name])
        except:
            print("Skipped SNP: %s" % snp_name)
            continue

        for upregulated_index, upregulated in upregulated_df.iterrows():
            gene_id = upregulated['ID']
            gene_upregulated_transitions = get_up_transitions(upregulated)

            upregulated_transitions = gene_upregulated_transitions - snp_gene_upregulated_transitions
            upregulated_transitions_str = ';'.join(upregulated_transitions)

            upregulated_transitions_str = multipleReplace(upregulated_transitions_str, snp_values)

            snp_dict[gene_id].update(
                allele_transition_counts(upregulated_transitions_str)
            )

        for gene_id in gene_ids:
            w.writerow(snp_dict[gene_id])
        print('-'*50)


def snp_ids_for_tfs():
    orf_manipulator = ORFManipulator()
    tfs = pandas.read_csv('../../../data/product/tf_snp_count.csv')
    tfs_with_snps = tfs[tfs['SNP_COUNT']>0]

    snps = []
    for tf_id in tfs_with_snps['ID']:
        snps = snps + list(orf_manipulator.orf_id_to_snp(tf_id).values)

    return snps


def create_tf_snp_count():
    orf_manipulator = ORFManipulator()
    snps = pandas.read_csv('../../../data/product/snp_orf.csv')
    tf_ids = list(pandas.read_csv('../../../data/product/annotated_true_tf_targets.csv').TF.values)

    # only SNPS in TFs
    snps = snps[snps.ID.apply(lambda x: x in tf_ids)]

    snps = snps.groupby('ID').count()
    snps['ID'] = snps.index
    snps['TF_NAME'] = snps.ID.apply(lambda x: orf_manipulator.orf_id_to_gene_name(x))
    snps = snps.rename(columns={'SNP_ID': 'SNP_COUNT'})
    snps = snps.dropna()

    snps = snps.sort_values(by=['SNP_COUNT'], ascending=False)
    snps.to_csv('../../../data/product/tf_snp_count.csv', index=False)
    return snps


# genes that are direct targets of SNP containing TFs
def create_tf_targets_file(): 
    tfs_with_snps = TFManipulator().tfs_with_snps
    tfs_all_targets = TargetManipulator().tf_targets

    df = tfs_all_targets.loc[tfs_with_snps['ID']]

    df = df.transpose()[df.apply(lambda x: x.sum()!=0, axis=0)]
    df = pandas.DataFrame(data=df.index, columns=['TARGET_ID'])

    df.to_csv('../../../data/product/targets_of_tfs_with_snps.csv', index=False)
    return df


# create an expression file with only the targets of TFs with SNPs and no missing values
def create_tf_expression_matrix():
    df = pandas.read_csv('../../../data/expression_matrix.txt', delimiter='\t', index_col='ID')
    genes = pandas.read_csv('../../../data/product/targets_of_tfs_with_snps.csv')

    # only targets on TFs
    df = df.iloc[genes['TARGET_ID']]
    # no missing values
    df = df.dropna()
    
    df = df.transpose()
    df.index.name = 'CONDITION'

    df.to_csv('../../../data/product/expression_matrix_tf_targets_only.csv')
    return df

