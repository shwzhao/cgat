#!/usr/bin/env python3

import argparse
from commands import gff2idmap, gffreorder, longestSeq, fa2agp, parse_orthofinder, taxonomy

def main():
    parser = argparse.ArgumentParser(prog='cgat', description='Comparative genomics analysis toolkit')
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand')

    gff2idmap.setup_parser(subparsers)
    gffreorder.setup_parser(subparsers)
    longestSeq.setup_parser(subparsers)
    fa2agp.setup_parser(subparsers)
    taxonomy.setup_parser(subparsers)
    parse_orthofinder.setup_parser(subparsers)


    args = parser.parse_args()

    if args.subcommand:
        if args.subcommand == 'gff2idmap':
            gff2idmap.run(args)
        elif args.subcommand == 'gffreorder':
            gffreorder.run(args)
        elif args.subcommand == 'longestSeq':
            longestSeq.run(args)
        elif args.subcommand == 'fa2agp':
            fa2agp.run(args)
        elif args.subcommand == 'taxonomy':
            taxonomy.run(args)
        elif args.subcommand == 'parse_orthofinder':
            parse_orthofinder.run(args)
        else:
            parser.print_help()
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
