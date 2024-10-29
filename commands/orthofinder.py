# commands/command2.py

def setup_parser(parser):
    orthofinder_parser = parser.add_parser('orthofinder', help='orthofinder help')
    # Add command 2 specific arguments
    orthofinder_parser.add_argument('arg2', help='Argument 2 for command 2')
    orthofinder_parser.add_argument('--option2', help='Option 2 for command 2')

def run(args):
    # Implement functionality for command 2
    print(f'Running command 2 with argument: {args.arg2}')
    if args.option2:
        print(f'Option 2 value: {args.option2}')
