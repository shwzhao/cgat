

# fasttree
# iq-tree
# raxml
# astral    https://github.com/smirarab/ASTRAL/ https://github.com/chaoszhang/ASTER https://doi.org/10.1093/bioinformatics/btu462

# paralog

def setup_parser(parser):
    parse_orthofinder_parser = parser.add_parser('parse_orthofinder', help='Parse orthofinder\'s results')
    # Add command 2 specific arguments
    parse_orthofinder_parser.add_argument('-f', '--orthofinder_result_file', required=True, help='Path to orthofinder result file')
    parse_orthofinder_parser.add_argument('-p', '--pattern', required=True, help='{align, retree, pangene}')
    parse_orthofinder_parser.add_argument('-o', '--option2', help='Option 2 for command 2')

    return parse_orthofinder_parser

def run(args):
    # Implement functionality for command 2
    print(f'Running command 2 with argument: {args.arg2}')
    if args.option2:
        print(f'Option 2 value: {args.option2}')
