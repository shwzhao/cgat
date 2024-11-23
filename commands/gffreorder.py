#!/usr/bin/env python3


def reverse_coordinates(start, end, chrom_length):
    """反转基因坐标"""
    new_start = chrom_length - end + 1
    new_end = chrom_length - start + 1
    return new_start, new_end

def update_gff(input_gff, rename_file, output_gff=None):
    # 读取 rename 文件，建立染色体名称映射表
    rename_dict = {}
    with open(rename_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            original, new_name, flip, chrom_length = parts[0], parts[1], parts[2] == "True", int(parts[3])
            rename_dict[original] = {'new_name': new_name, 'flip': flip, 'length': chrom_length}
    
    # 打开输出文件或使用标准输出
    if output_gff is None:
        import sys
        out_f = sys.stdout  # 使用标准输出
    else:
        out_f = open(output_gff, 'w')  # 打开指定的输出文件

    try:
        # 读取并处理 GFF 文件
        with open(input_gff, 'r') as in_f:
            for line in in_f:
                if line.startswith("#"):
                    # 保留注释行
                    out_f.write(line)
                    continue
                
                fields = line.strip().split('\t')
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                
                # 判断是否需要改变染色体名称和方向
                if chrom in rename_dict:
                    new_name = rename_dict[chrom]['new_name']
                    flip = rename_dict[chrom]['flip']
                    chrom_length = rename_dict[chrom]['length']
                    
                    # 如果需要翻转方向，更新坐标和方向
                    if flip:
                        new_start, new_end = reverse_coordinates(start, end, chrom_length)
                        start, end = new_start, new_end
                        strand = "-" if strand == "+" else "+"
                    
                    # 更新染色体名称
                    chrom = new_name
                
                # 更新并写入 GFF 文件
                fields[0] = chrom
                fields[3] = str(start)
                fields[4] = str(end)
                fields[6] = strand
                out_f.write('\t'.join(fields) + '\n')
    finally:
        # 确保文件在不为 stdout 时正确关闭
        if output_gff is not None:
            out_f.close()

def setup_parser(parser):
    gffreorder_parser = parser.add_parser('gffreorder', help='Reorder gff file according to ref')

    gffreorder_parser.add_argument('-g', '--input_gff', required=True, help='Path to gff file')
    gffreorder_parser.add_argument('-o', '--output_gff', help='Path to the output file.')
    gffreorder_parser.add_argument('-r', '--rename_file', help='rename file, #q_id r_id ')
    
    return gffreorder_parser

def run(args):
    update_gff(args.input_gff, args.rename_file, args.output_gff)
