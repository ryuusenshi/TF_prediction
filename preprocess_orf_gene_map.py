import pandas

HEADER_LINE_COUNT = 27

df = pandas.read_csv('data/GPL118.annot', delimiter='\t', header=HEADER_LINE_COUNT)
if df.columns.values[0] != 'ID':
    raise ValueError('Expected ID as the first column of the microarray ORF ID annotation file: "data/GPL118.annot". Did the number of header lines change?')

# remove unnecessary columns
df = pandas.DataFrame(data=df, columns=['ID', 'Gene title', 'Gene symbol', 'Platform_ORF'])
# remove probes that are not testing for ORFs
df = df.loc[df.Platform_ORF.dropna().index]
df.to_csv('data/id_to_gene_name.csv', index=False)
