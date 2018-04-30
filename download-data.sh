# make data directory
mkdir data
# make data intermediate product directory
mkdir data/product

# downloads the Xie et al. preprocessed data
wget https://github.com/shilab/MLP-SAE/archive/master.zip
unzip master.zip
mv MLP-SAE-master/data/* ./data
rm -rf MLP-SAE-master
rm master.zip

# downloads the known TF and known targets interaction matrix data
wget http://www.yeastract.com/download/RegulationMatrix_Documented_2013927.csv.gz
gunzip RegulationMatrix_Documented_2013927.csv.gz
mv RegulationMatrix_Documented_2013927.csv data/RegulationMatrix.csv

# downloads the microarray ORF ID annotation file
wget ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL118/annot/GPL118.annot.gz
gunzip GPL118.annot.gz
mv GPL118.annot data/

# transforms the microarray ORD ID annotation file from SOFT format to csv
python preprocess_orf_gene_map.py
rm data/GPL118.annot
