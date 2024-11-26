#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def parse_fasta(fasta):
    '''
    https://github.com/zengxiaofei/HapHiC/blob/main/utils/mock_agp_file.py
    # Author: Xiaofei Zeng
    # Email: xiaofei_zeng@whu.edu.cn
    # Created Time: 2021-10-22 11:20
    '''
    import collections
    from commands.read_data import open_file

    len_dict = collections.OrderedDict()
    with open_file(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                len_dict[ID] = 0
            else:
                len_dict[ID] += len(line.strip())
    return len_dict


def mock_agp(len_dict, output_file=None):
    from commands.read_data import write_output

    output_lines = []
    for ID, length in len_dict.items():
        line = '{0}\t1\t{1}\t1\tW\t{0}\t1\t{1}\t+\n'.format(ID, length)
        output_lines.append(line)

    write_output(output_lines, output_file)


def generate_agp(input_file, output_file, gap_size, contig_prefix, gap_identifer="N"):
    """
    将基因组组装的 FASTA 文件转换为 AGP 文件。
    :param fasta_file: 输入的 FASTA 文件路径。
    :param output_file: 输出的 AGP 文件路径。
    :param gap_size: 定义 gap 的最小长度（默认为连续 N 的数量）。
    :param contig_prefix: contig 名称的前缀。
    """
    from Bio import SeqIO
    from commands.read_data import open_file, write_output

    output_lines = []
    gap_identifer=gap_identifer.upper()
    with open_file(input_file) as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            chrom_name = record.id  # 染色体名称
            sequence = record.seq.upper()  # 将序列转换为大写
            start = 1
            contig_count = 1

            # 遍历染色体序列，识别 gaps 和 contigs
            i = 0
            while i < len(sequence):
                if sequence[i] == gap_identifer:  # 进入 gap 区域
                    gap_start = i
                    while i < len(sequence) and sequence[i] == gap_identifer:
                        i += 1
                    gap_end = i
                    gap_length = gap_end - gap_start
                    if gap_length >= gap_size:  # 符合 gap 条件
                        # 输出 gap 信息到 AGP 文件
                        line = f"{chrom_name}\t{start}\t{start + gap_length - 1}\t{contig_count}\tU\t{gap_length}\tscaffold\tyes\tproximity_ligation\n"
                        start += gap_length
                        # contig_count += 1
                        output_lines.append(line)
                    continue
                else:  # 进入 contig 区域
                    contig_start = i
                    while i < len(sequence) and sequence[i] != gap_identifer:
                        i += 1
                    contig_end = i
                    contig_length = contig_end - contig_start
                    # 输出 contig 信息到 AGP 文件
                    line = f"{chrom_name}\t{start}\t{start + contig_length - 1}\t{contig_count}\tW\t{contig_prefix}{contig_count}\t1\t{contig_length}\t+\n"
                    start += contig_length
                    contig_count += 1
                    output_lines.append(line)

    write_output(output_lines, output_file)

def setup_parser(parser):
    fa2agp_parser = parser.add_parser('fa2agp', help='Convert genome assembly FASTA format to AGP file')
    # Add command specific arguments
    fa2agp_parser.add_argument('-i', '--input_file', required=True, help='Path to the input fasta file.')
    fa2agp_parser.add_argument('-o', '--output_file', help='Path to the output agp file. [sdout]')
    fa2agp_parser.add_argument('-m', '--min_gap_size', type=int, default=10, help='Minimum length of gap to split contigs. [10]')
    fa2agp_parser.add_argument('-n', '--gap_identifer', default="N", help='Gap identifer (ignore case). [N]')
    fa2agp_parser.add_argument('-p', '--contig_prefix', type=str, default="contig_", help='Prefix for new contig names. [contig_]')

    return fa2agp_parser

def run(args):
    if args.min_gap_size == 0:
        len_dict = parse_fasta(args.input_file)
        mock_agp(len_dict)
    else:
        generate_agp(args.input_file, args.output_file, args.min_gap_size, args.contig_prefix, args.gap_identifer)
