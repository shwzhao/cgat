#!/usr/bin/env python3
'''
除了基本的信息, 还可以添加其他信息, 比如gene行的Description
还可以提取其他行的ID和parent, 比如有的gff没有mRNA行, 可以提取Transcript行的ID和Parent
gene::Description;Transcript::ID;Transcript::Parent
'''

def parse_extra_columns(extra_columns_str):
    extra_columns = extra_columns_str.split(';')
    column_mapping = {}
    for column in extra_columns:
        key, value = column.split('::')
        if key not in column_mapping:
             column_mapping[key] = []
        column_mapping[key].append(value)
    return extra_columns, column_mapping


def parse_gff(gff_file, mRNA_Type = 'mRNA', extra_columns = ''):
    gene_id_mapping = {}
    mrna_id_mapping = {}
    cds_id_mapping = {}

    if not extra_columns:
        extra_columns = ''
        extra_columns_mapping = {}
    else:
        extra_columns, extra_columns_mapping = parse_extra_columns(extra_columns)

    from commands.read_data import open_file
    with open_file(gff_file) as f:
        for line in f:
            fields = line.strip().strip('\r;').split('\t')
            if line.startswith('#'):
                continue
            if len(fields) != 9:
                continue
            SeqID, Source, Type, Start, End, Score, Strand, Phase, Attributes = fields

            attr_dict = {}
            Attributes = Attributes.split(';')
            for attr in Attributes:
                key_value = attr.strip().split('=')
                if len(key_value) == 2:  # 忽略不完整的键值对
                    key, value = key_value
                    attr_dict[key] = value

            gene_extra_value_dict = {}
            rna_extra_value_dict = {}

            if Type == 'gene':
                gene_id = attr_dict.get('ID', '')
                gene_name = attr_dict.get('Name', gene_id)
                gene_id_mapping[gene_id] = {}
                gene_id_mapping[gene_id]['gene_name'] = gene_name
                if Type in extra_columns_mapping.keys():
                    extra_values = [attr_dict.get(extra_attr, None) for extra_attr in extra_columns_mapping[Type]]
                    gene_extra_value_pairs = zip(['Extra_' + Type + '_' + value for value in extra_columns_mapping[Type]], extra_values)
                    gene_extra_value_dict = {key: value for key, value in gene_extra_value_pairs}
                    gene_id_mapping[gene_id].update(gene_extra_value_dict)

            elif Type == mRNA_Type:
                rna_id = attr_dict.get('ID', '')
                gene_id = attr_dict.get('Parent', rna_id)
                transcript_name = attr_dict.get('Name', rna_id)
                mrna_id_mapping[rna_id] = {}
                mrna_id_mapping[rna_id]['gene_id'] = gene_id
                mrna_id_mapping[rna_id]['transcript_id'] = rna_id
                mrna_id_mapping[rna_id]['transcript_name'] = transcript_name
                mrna_id_mapping[rna_id]['SeqID'] = SeqID
                mrna_id_mapping[rna_id]['Start'] = Start
                mrna_id_mapping[rna_id]['End'] = End
                mrna_id_mapping[rna_id]['Strand'] = Strand

                if Type in extra_columns_mapping.keys():
                    extra_values = [attr_dict.get(extra_attr, None) for extra_attr in extra_columns_mapping[Type]]
                    rna_extra_value_pairs = zip(['Extra_' + Type + '_' + value for value in extra_columns_mapping[Type]], extra_values)
                    rna_extra_value_dict = {key: value for key, value in rna_extra_value_pairs}
                    mrna_id_mapping[rna_id].update(rna_extra_value_dict)

            elif Type == "CDS":
                cds_id = attr_dict.get('ID', 'Parent')
                cds_id_mapping[attr_dict['Parent']] = cds_id

            else:
                continue
    return gene_id_mapping, mrna_id_mapping, cds_id_mapping, extra_columns

def write_idmapping_file(gene_id_mapping, mrna_id_mapping, cds_id_mapping, extra_columns, output_file):
    with open(output_file, 'w') as f:
        # f.write("#GeneID\tGeneName\n")
        keys_to_output_first = ['gene_id', 'gene_name', 'transcript_id', 'transcript_name', "cds_id", 'SeqID', 'Start', 'End', 'Strand']
        if extra_columns:
            header_line = '#' + '\t'.join(keys_to_output_first + extra_columns)
        else:
            header_line = '#' + '\t'.join(keys_to_output_first)
        f.write(f"{header_line}\n")
        for rna_id, rna_attr_dict in mrna_id_mapping.items():
            gene_id = rna_attr_dict['gene_id']
            try:
                rna_attr_dict.update(gene_id_mapping[gene_id])
            except KeyError:
                print(f"Warning: {gene_id} not found in gene_id_mapping. Skipping update.")
            rna_attr_dict["cds_id"] = cds_id_mapping.get(rna_id, None)

            output_line = '\t'.join(str(rna_attr_dict.get(key, None)) for key in keys_to_output_first)
            output_line += '\t' + '\t'.join(str(value) for key, value in reversed(rna_attr_dict.items()) if key not in keys_to_output_first)
            # print(rna_attr_dict)
            
            f.write(f"{output_line}\n")

def setup_parser(parser):
    idmap_parser = parser.add_parser('gff2idmap', help='Convert GFF file to ID_MAP format')
    
    idmap_parser.add_argument('-g', '--gff_file', required=True, help='Path to gff file')
    idmap_parser.add_argument('-o', '--output_file', default='id_mapping.txt', help='Path to the output file. [id_mapping.txt]')
    idmap_parser.add_argument('-t', '--trans_mRNA_info_to', default='mRNA', help='Transcript or mRNA. [mRNA]')
    idmap_parser.add_argument('-e', '--extra_info', help='Extra information that you need, for example: -e "mRNA::Dbxref;gene::gbkey". [NULL]')
    idmap_parser.add_argument('--only_coding_gene', action='store_true', help='only map pep coding gene ID')

    return idmap_parser
    '''
    如果-e有除了gene和mRNA的情况, 需要同时改变-t参数, 如果改变了-t参数, 就不能再得到mRNA的其他信息
    '''

def run(args):
    # Implement functionality for command 1
    gene_id_mapping, mrna_id_mapping, cds_id_mapping, extra_columns = parse_gff(args.gff_file, mRNA_Type = args.trans_mRNA_info_to, extra_columns=args.extra_info)
    write_idmapping_file(gene_id_mapping, mrna_id_mapping, cds_id_mapping, extra_columns, args.output_file)
