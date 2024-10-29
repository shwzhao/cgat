
# ggtree


# https://github.com/tanghaibao/jcvi
# https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version)





def setup_parser(parser):
    tree_parser = parser.add_parser('tree', help='get longest transcript help')
    # Add command specific arguments
    tree_parser.add_argument('-i', '--idmapping_file', required=True, help='Path to the gene-transcript mapping file')
    tree_parser.add_argument('-s', '--transcript_file', required=True, help='Path to the transcript sequences file')
    tree_parser.add_argument('-o', '--output_file', default='output.fa', help='Path to the output file. | Default: output.fa')
    tree_parser.add_argument('-l', '--length_file', help='Path to the output transcript length file')
    tree_parser.add_argument('-n', '--number',  type=int, default=2, help='Which column you want to map. | Default: 2')
    tree_parser.add_argument('-d', '--not_change_name', action='store_true', help='Do not change transcript name to gene name.')

    return tree_parser

def run(args):
    # 读取基因名和转录本名的对应关系
    print("hello world!")




# dotplot图
# 基因组共线性图


# 简单的统计的图




