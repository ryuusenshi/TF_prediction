import csv
import numpy
import re
import pandas

from collections import defaultdict


# log threshold for up/down regulation
log_threshold = 1


class ORFManipulator(object):
    _orf_to_snp_df = None

    def __init__(self):
        if self._orf_to_snp_df is None:
            self._orf_to_snp_df = pandas.read_csv('snp_orf.csv')
            self._orf_info = pandas.read_csv('id_to_gene_name.csv')
 
    def orf_id_to_snp(self, orf_id):
        df = self._orf_to_snp_df.loc[self._orf_to_snp_df['ID'] == int(orf_id)]

        return df['SNP_ID']

    def snp_to_orf_id(self, snp_id):
        df = self._orf_to_snp_df.loc[self._orf_to_snp_df['SNP_ID'] == int(snp_id)]

        return df['ID'].values[0]

    def orf_id_to_orf_name(self, orf_id):
        df = self._orf_info.loc[self._orf_info['ID'] == int(orf_id)]

        return df['ORF_MA'].values[0]

    def orf_name_to_orf_id(self, orf_name):
        df = self._orf_info.loc[self._orf_info['ORF_MA'] == orf_name]

        return df['ID'].values[0]

    def orf_id_to_gene_name(self, orf_id):
        df = self._orf_info.loc[self._orf_info['ID'] == int(orf_id)]

        return df['GENE'].values[0]



orf_manipulator = ORFManipulator()


def create_snp_orf_file():
    # inputs
    orf_df = pandas.read_csv('expression_location', delimiter='\t')
    snp_df = pandas.read_csv('snp_position.txt', delimiter='\t')

    # output
    f = open('snp_orf.csv', 'w')
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

def create_snp_target_file():
    # only TFs that have at least one SNP
    df = pandas.read_csv('tf_snp_count.csv', delimiter=',')
    df = df.loc[df['SNP_COUNT']>0]

    for index, tf in df.iterrows():
        orf_id = tf['ID']


def create_upregulated_file_v2():
    # input file
    exp_df = pandas.read_csv('expression_matrix.txt', delimiter='\t')

    # output file
    # output file
    out_filename = 'upregulated_targets_v2.txt'
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

def create_upregulated_file():
    """
        inputs: expression_matrix.txt
        ouputs: upregulated_targets.txt
    """
    # input file
    in_filename = 'expression_matrix.txt'
    r = csv.DictReader(open(in_filename), delimiter='\t')
    d = numpy.genfromtxt(open(in_filename, "rb"), delimiter="\t")

    # output file
    out_filename = 'upregulated_targets.txt'
    f = open(out_filename, 'w')
    w = csv.writer(f)
    # write header to output file
    w.writerow(['ID', 'sample1', 'sample2'])
    

    for index, row in enumerate(d[1:]):
        ID = '%d' % row[0]
        n = len(row)-1
        total_pairs = n*(n-1)
        row = row[1:]
        for i in range(n):
            for j in range(i+1,n):
                if row[i]-row[j]>=log_threshold:
                    cond1 = r.fieldnames[i+1]
                    cond2 = r.fieldnames[j+1]
                if row[j]-row[i]>=log_threshold:
                    cond1 = r.fieldnames[j+1]
                    cond2 = r.fieldnames[i+1]
                if abs(row[i]-row[j])>=log_threshold:
                    w.writerow([ID, cond1, cond2])
    f.close()

def create_snp_to_gene_file_v2():
    # -- is shorthand for 00, 11, and 22 (no change)
    # ^^ is shorhand for 02, 01, 12 (minor to major)
    # vv is shorthand for 20, 10, 21 (major to minor)
    # IMPORTANT: SNP IDS are counted from 0
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

    filename = 'snp_to_gene_v2.txt'
    f = open(filename, 'w')
    fieldnames = ['SNP_ID', 'GENE_ID', '0--2', '0--1', '1--2', '2--0', '1--0', '2--1', '0--0', '1--1', '2--2'] 
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()


    # only SNPs for 5513 TF
    snps = orf_manipulator.orf_id_to_snp(5513)
    snp_df = pandas.read_csv('matrix_genotypes.txt', delimiter='\t')

    
    snp_df = snp_df.loc[snps]
    conditions = snp_df.columns.values[1:]
    # add SNP_ID to the snp_df dataframe
    snps.index = snp_df.index
    snp_df = snp_df.join(snps)


    upregulated_df = pandas.read_csv('upregulated_targets_v2.txt')
    upregulated_df = upregulated_df.dropna()
    for _, snp in snp_df.iterrows(): 
        snp_name = snp['SNP_ID']
        snp_values = snp.to_dict()
        print("Working on SNP: %d" % snp_name)
        # TODO: remove hardcode
        snp_gene_id = 5513

        snp_dict = {}
        gene_ids = set(upregulated_df['ID']) 
        for gene_id in gene_ids:
            snp_dict[gene_id] = defaultdict(int)
            snp_dict[gene_id].update({
                'SNP_ID': snp_name,
                'GENE_ID': gene_id,
            })
       
        snp_gene_upregulated_transitions = get_up_transitions(upregulated_df.loc[snp_name])
        for upregulated_index, upregulated in upregulated_df.iterrows():
            gene_id = upregulated['ID']
            gene_upregulated_transitions = get_up_transitions(upregulated)

            upregulated_transitions = gene_upregulated_transitions - snp_gene_upregulated_transitions
            upregulated_transitions_str = ';'.join(upregulated_transitions)

            upregulated_transitions_str = multipleReplace(upregulated_transitions_str, snp_values)

            snp_dict[gene_id].update(
                allele_transition_counts(upregulated_transitions_str)
            )

            if not upregulated_index % 200000 and upregulated_index > 0:
                print('  ---- processed %d upregulated pairs' % upregulated_index)
        for gene_id in gene_ids:
            w.writerow(snp_dict[gene_id])
        print('-'*50)

def create_snp_to_gene_file():
    # -- is shorthand for 00, 11, and 22 (no change)
    # ^^ is shorhand for 02, 01, 12 (minor to major)
    # vv is shorthand for 20, 10, 21 (major to minor)
    # IMPORTANT: SNP IDS are counted from 0
    filename = 'snp_to_gene.txt'
    f = open(filename, 'w')
    fieldnames = ['SNP_ID', 'GENE_ID', '--', '02', '01', '12', '20', '10', '21', '00', '11', '22'] 
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()


    snp_df = pandas.read_csv('matrix_genotypes.txt', delimiter='\t')
    # snp_df = snp_df.iloc[68:69]
    upregulated_df = pandas.read_csv('upregulated.txt')
    for row_index in range(len(snp_df)):
        print("Working on SNP: %d" % row_index)
        row = snp_df.iloc[row_index]
        snp_name = row_index
        snp_dict = {}
        gene_ids = set(upregulated_df['ID']) 
        for gene_id in gene_ids:
            snp_dict[gene_id] = defaultdict(int)
            snp_dict[gene_id].update({
                'SNP_ID': snp_name,
                'GENE_ID': gene_id,
            })
            
        for upregulated_index in range(len(upregulated_df)):
            upregulated_row = upregulated_df.iloc[upregulated_index]
            gene_id = upregulated_row['ID']
            upsample_name = upregulated_row['sample1']
            downsample_name = upregulated_row['sample2']

            upsample_snp = row[upsample_name]
            downsample_snp = row[downsample_name]

            transition_id = '%s%s' % (downsample_snp, upsample_snp)
            snp_dict[gene_id][transition_id] += 1

            if downsample_snp == upsample_snp:
                snp_dict[gene_id]['--'] += 1
#
#            if upsample_snp == downsample_snp:
#                snp_dict[gene_id]['SAMESNP'] += 1
#            else:
#                snp_dict[gene_id]['UPSNP%s' % upsample_snp] += 1
#                snp_dict[gene_id]['DOWNSNP%s' % downsample_snp] += 1
            if not upregulated_index % 200000 and upregulated_index > 0:
                print('  ---- processed %d upregulated pairs' % upregulated_index)
        for gene_id in gene_ids:
            w.writerow(snp_dict[gene_id])
        print('-'*50)
