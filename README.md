# CGAT

## 1. Description
CGAT (Comparative Genomics Analysis Toolkit)

## 2. Installation
```
git clone https://github.com/shwzhao/cgat
```

## 3. Usage
```

```

### gff2idmap
```
cgat gff2idmap -h
usage: main.py gff2idmap [-h] -g GFF_FILE [-o OUTPUT_FILE] [-t TRANS_MRNA_INFO_TO] [-e EXTRA_INFO] [--only_coding_gene]

options:
  -h, --help            show this help message and exit
  -g GFF_FILE, --gff_file GFF_FILE
                        Path to gff file
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Path to the output file. [id_mapping.txt]
  -t TRANS_MRNA_INFO_TO, --trans_mRNA_info_to TRANS_MRNA_INFO_TO
                        Transcript or mRNA. [mRNA]
  -e EXTRA_INFO, --extra_info EXTRA_INFO
                        Extra information that you need, for example: -e "mRNA::Dbxref;gene::gbkey". [NULL]
  --only_coding_gene    only map pep coding gene ID
```



```
cgat gff2idmap -g example/GCF_000002775.5_P.trichocarpa_v4.1_genomic.gff.gz

head -5 id_mapping.txt
## #gene_id        gene_name       transcript_id   transcript_name cds_id  SeqID   Start   End     Strand
## gene-LOC7483226 LOC7483226      rna-XM_002299051.3      XM_002299051.3  cds-XP_002299087.1      NC_037285.2     40829   43223   -
## gene-LOC7457680 LOC7457680      rna-XM_002297573.4      XM_002297573.4  cds-XP_002297609.3      NC_037285.2     50624   52700   +
## gene-LOC7483220 LOC7483220      rna-XM_006368463.3      XM_006368463.3  cds-XP_006368525.3      NC_037285.2     60500   65314   -
## gene-LOC7457681 LOC7457681      rna-XM_024603583.2      XM_024603583.2  cds-XP_024459351.2      NC_037285.2     65961   72394   -
```

You can also output extra information you want:
```
cgat gff2idmap -g example/GCF_000002775.5_P.trichocarpa_v4.1_genomic.gff.gz -e "gene::Dbxref;mRNA::experiment"
head -5 id_mapping.txt 
## #gene_id        gene_name       transcript_id   transcript_name cds_id  SeqID   Start   End     Strand  gene::Dbxref    mRNA::experiment
## gene-LOC7483226 LOC7483226      rna-XM_002299051.3      XM_002299051.3  cds-XP_002299087.1      NC_037285.2     40829   43223   -       GeneID:7483226  COORDINATES: polyA evidence [ECO:0006239]
## gene-LOC7457680 LOC7457680      rna-XM_002297573.4      XM_002297573.4  cds-XP_002297609.3      NC_037285.2     50624   52700   +       GeneID:7457680  COORDINATES: polyA evidence [ECO:0006239]
## gene-LOC7483220 LOC7483220      rna-XM_006368463.3      XM_006368463.3  cds-XP_006368525.3      NC_037285.2     60500   65314   -       GeneID:7483220  COORDINATES: polyA evidence [ECO:0006239]
## gene-LOC7457681 LOC7457681      rna-XM_024603583.2      XM_024603583.2  cds-XP_024459351.2      NC_037285.2     65961   72394   -       GeneID:7457681  COORDINATES: polyA evidence [ECO:0006239]
```

### longestSeq

Using the id_mapping file to generate gene sequence file that there are only longest isoform. You can also change the gene name.
```
cgat longestSeq -h
usage: main.py longestSeq [-h] -i IDMAPPING_FILE -s TRANSCRIPT_FILE [-o OUTPUT_FILE] [-l LENGTH_FILE]
                          [-n NUMBER] [-d]

options:
  -h, --help            show this help message and exit
  -i IDMAPPING_FILE, --idmapping_file IDMAPPING_FILE
                        Path to the gene-transcript mapping file
  -s TRANSCRIPT_FILE, --transcript_file TRANSCRIPT_FILE
                        Path to the transcript sequences file
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Path to the output file. | Default: output.fa
  -l LENGTH_FILE, --length_file LENGTH_FILE
                        Path to the output transcript length file
  -n NUMBER, --number NUMBER
                        Which column you want to map. | Default: 2
  -d, --not_change_name
                        Do not change transcript name to gene name.
```

```
zcat example/GCF_000002775.5_P.trichocarpa_v4.1_protein.faa.gz | grep ">" | head -3
## >NP_001391697.1 caffeoyl-CoA O-methyltransferase 1 [Populus trichocarpa]
## >NP_001391698.1 caffeoyl-CoA O-methyltransferase 2 [Populus trichocarpa]
## >XP_002297609.3 L-ascorbate oxidase homolog [Populus trichocarpa]

# There no matched names in id_mapping.txt, you can change the id_mapping.txt column or fasta sequence name to fit.

awk 'BEGIN{OFS="\t"}{gsub(/cds-/, "", $5);print}' id_mapping.txt \
  | cgat longestSeq \
  -i - \
  -s example/GCF_000002775.5_P.trichocarpa_v4.1_protein.faa.gz \
  -o pep.faa \
  -n 5

grep ">" pep.faa | head -3
## >gene-LOC7483226
## >gene-LOC7457680
## >gene-LOC7483220
```


