import csv

from collections import defaultdict


def tf_to_id_table():
    f = open('tf_name_to_id.csv', 'w')
    f_missing = open('tf_missing_ids.csv', 'w')

    reader = csv.DictReader(open('id_to_gene_name.csv'))
    writer = csv.DictWriter(f, fieldnames=['TF_NAME', 'ID'])
    writer.writeheader()

    known_tfs = open('known_tfs.csv').read().strip().split('\n')

    # gene name to microarray id mapping
    gene_dict = dict()
    orf_dict = dict()
    for row in reader:
        gene_dict[row['GENE'].lower()] = row['ID']
        orf_dict[row['ORF']] = row['ID']


    for known_tf in known_tfs:
        # true name is lowercase and might exclude the trailing letter 'p'
        known_tf_true_name = known_tf.lower()
        try:
            gene_id = gene_dict[known_tf_true_name]
        except KeyError:
            known_tf_true_name = known_tf_true_name[:-1]
            try:
                gene_id = gene_dict[known_tf_true_name]
            except KeyError:
                try:
                    gene_id = orf_dict[known_tf]
                except KeyError:
                    print("Transcription factor: %s not found in the microarray" % known_tf)
                    f_missing.write('%s\n' % known_tf)
                    continue

        writer.writerow({
            'ID': gene_id,
            'TF_NAME': known_tf,
        })


def tf_snp_count():
    reader = csv.reader(open('ID_Chr_SNP.csv'))
    gene_dict = defaultdict(int)
    for row in reader:
        gene_id = row[0]
        gene_dict[gene_id] += 1

    f = open('tf_snp_count.csv', 'w')
    reader = csv.DictReader(open('tf_name_to_id.csv'))
    writer = csv.DictWriter(f, fieldnames=['ID', 'TF_NAME', 'SNP_COUNT'])
    writer.writeheader()

    output_rows = []
    for row in reader:
        gene_id = row['ID']
        snp_count = gene_dict[gene_id]

        output_rows.append({
            'ID': row['ID'],
            'TF_NAME': row['TF_NAME'],
            'SNP_COUNT': snp_count,
        })
    # sort by SNP count
    output_rows = sorted(output_rows, key=lambda d: d['SNP_COUNT'], reverse=True)

    # write all rows
    for output_row in output_rows:
        writer.writerow(output_row)
    f.close()


def remove_unknowns():
    f = open('id_to_gene_name_no_unknowns.csv', 'w')
    reader = csv.DictReader(open('id_to_gene_name.csv'))
    writer = csv.DictWriter(f, fieldnames=reader.fieldnames)
    writer.writeheader()

    for row in reader:
        if row['ORF'].lower() != 'unknown':
            writer.writerow(row)

    f.close()
