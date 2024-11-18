from Bio import SeqIO

def read_gene_transcript_mapping(mapping_file, map_column = 2):
    """
    从基因名和转录本名的对应关系文件中读取基因名和转录本名的对应关系。
    """
    gene_transcripts = {}
    from commands.read_data import open_file
    with open_file(mapping_file) as f:
    # with open(mapping_file, 'r') as f:
        for line in f:
            name_list = line.strip().split('\t')
            gene_name, transcript_name = name_list[0], name_list[map_column-1]
            if gene_name not in gene_transcripts:
                gene_transcripts[gene_name] = []
            gene_transcripts[gene_name].append(transcript_name)
    return gene_transcripts

def read_transcript_sequences(transcript_file):
    """
    从转录本序列文件中读取转录本序列，并返回转录本名和序列的对应关系。
    """
    transcript_sequences = {}
    from commands.read_data import open_file
    with open_file(transcript_file) as f:
    # with open(transcript_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            transcript_name = record.id
            sequence = str(record.seq)
            transcript_sequences[transcript_name] = sequence
    return transcript_sequences

def find_longest_transcripts(gene_transcripts, transcript_sequences, change_name = True, output_file = None):
    """
    根据每个基因名对应的转录本长度信息，找到每个基因名对应的最长转录本序列。
    注意: 并不一定每个ID都对应了序列
    """
    longest_transcripts = {}
    if output_file:
        with open(output_file, 'w') as f:
            for gene_name, transcripts in gene_transcripts.items():
                longest_transcript = max(transcripts, key=lambda x: len(transcript_sequences[x]))
                longest_transcripts[gene_name] = transcript_sequences[longest_transcript]
                if change_name:
                    longest_transcripts[longest_transcript] = longest_transcripts.pop(gene_name)
                for transcript in transcripts:
                    f.write(f'{gene_name}\t{transcript}\t{len(transcript_sequences[transcript])}\n')
    else:
        for gene_name, transcripts in gene_transcripts.items():
            try:
                longest_transcript = max(transcripts, key=lambda x: len(transcript_sequences[x]))
                longest_transcripts[gene_name] = transcript_sequences[longest_transcript]
                if change_name:
                    longest_transcripts[longest_transcript] = longest_transcripts.pop(gene_name)
            except KeyError:
                continue
    return longest_transcripts

def write_longest_transcript_sequences(longest_transcripts, output_file):
    """
    将每个基因名对应的最长转录本序列写入到输出文件中。
    """
    with open(output_file, 'w') as f:
        for gene_name, sequence in longest_transcripts.items():
            f.write(f'>{gene_name}\n{sequence}\n')

def setup_parser(parser):
    longestSeq_parser = parser.add_parser('longestSeq', help='get longest transcript help')
    # Add command specific arguments
    longestSeq_parser.add_argument('-i', '--idmapping_file', required=True, help='Path to the gene-transcript mapping file')
    longestSeq_parser.add_argument('-s', '--transcript_file', required=True, help='Path to the transcript sequences file')
    longestSeq_parser.add_argument('-o', '--output_file', default='output.fa', help='Path to the output file. | Default: output.fa')
    longestSeq_parser.add_argument('-l', '--length_file', help='Path to the output transcript length file')
    longestSeq_parser.add_argument('-n', '--number',  type=int, default=2, help='Which column you want to map. | Default: 2')
    longestSeq_parser.add_argument('-d', '--not_change_name', action='store_true', help='Do not change transcript name to gene name.')

    return longestSeq_parser

def run(args):
    # 读取基因名和转录本名的对应关系
    gene_transcripts = read_gene_transcript_mapping(args.idmapping_file, map_column=args.number)
    # 读取转录本序列
    transcript_sequences = read_transcript_sequences(args.transcript_file)
    # 找到每个基因名对应的最长转录本
    longest_transcripts = find_longest_transcripts(gene_transcripts, transcript_sequences, change_name=args.not_change_name, output_file = args.length_file)
    # 将最长转录本序列写入到输出文件中
    write_longest_transcript_sequences(longest_transcripts, args.output_file)

