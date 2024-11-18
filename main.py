#!/usr/bin/env python3

import argparse
from commands import gff2idmap, longestSeq, orthofinder, taxonomy, renamegff

def main():
    parser = argparse.ArgumentParser(description='Comparative genomics analysis toolkit')
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand')

    gff2idmap.setup_parser(subparsers)
    longestSeq.setup_parser(subparsers)
    taxonomy.setup_parser(subparsers)
    orthofinder.setup_parser(subparsers)
    renamegff.setup_parser(subparsers)

    args = parser.parse_args()

    if args.subcommand:
        if args.subcommand == 'gff2idmap':
            gff2idmap.run(args)
        elif args.subcommand == 'longestSeq':
            longestSeq.run(args)
        elif args.subcommand == 'taxonomy':
            taxonomy.run(args)
        elif args.subcommand == 'orthofinder':
            orthofinder.run(args)
        elif args.subcommand == 'renamegff':
            renamegff.run(args)
        else:
            parser.print_help()
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
