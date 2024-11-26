import os
import pandas as pd
# import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from itertools import combinations
# from collections import defaultdict

class OrthoFinderResults:
    def __init__(self, result_dir):
        self.result_dir = result_dir
        self.gene_count_file = os.path.join(result_dir, "Orthogroups/Orthogroups.GeneCount.tsv")
        self.orthogroups_file = os.path.join(result_dir, "Orthogroups/Orthogroups.tsv")

    def load_gene_counts(self):
        """Load gene count data."""
        return pd.read_csv(self.gene_count_file, sep="\t", index_col=0).iloc[:, 0:-1]

    def load_orthogroups(self):
        """Load orthogroups data."""
        return pd.read_csv(self.orthogroups_file, sep="\t", index_col=0)

    def classify_gene_families(self, threshold_pct=0.9, threshold_count=None):
        """
        Classify gene families into core, softcore, dispensable, and private.
        Either threshold_pct or threshold_count must be set.
        """
        gene_counts = self.load_gene_counts()
        num_species = len(gene_counts.columns)
        # non_zero_counts = (gene_counts > 0).sum(axis=1)
        non_zero_counts = (gene_counts > 0).astype(int)

        if threshold_count:
            threshold = threshold_count
        elif threshold_pct:
            import math
            threshold = math.ceil(num_species * threshold_pct) # threshold_pct=0.9
        else:
            raise ValueError("Either threshold_pct or threshold_count must be set.")

        print("softcore threshold is: {}.".format(threshold))

        def classify(count):
            total_columns = gene_counts.shape[1]
            if count == total_columns:
                return "core"
            elif count >= threshold:
                return "softcore"
            elif count > 1:
                return "dispensable"
            else:
                return "private"

        non_zero_counts['category'] = non_zero_counts.sum(axis=1).map(classify)
        return non_zero_counts

    def extract_genes(self, non_zero_counts):
        """Extract genes by classification from Orthogroups.tsv."""
        orthogroups = self.load_orthogroups()
        category_df = non_zero_counts['category'].reset_index()

        long_df = orthogroups.reset_index().melt(id_vars=['Orthogroup'], var_name='Species', value_name='Genes')
        long_df = long_df.dropna(subset=['Genes'])
        merged_df = pd.merge(long_df, category_df, on='Orthogroup', how='left')
        merged_df = merged_df[['Orthogroup', 'category', 'Species', 'Genes']]

        return merged_df

    def core_pan_analysis(self, max_combinations=1000):
        """
        Perform core and pan gene analysis.
        """
        gene_counts = self.load_gene_counts()
        samples = gene_counts.columns
        pan_core_results = []

        for i in range(1, len(samples) + 1):
            comb = list(combinations(samples, i))
            if len(comb) > max_combinations:
                comb = comb[:max_combinations]
            for subset in comb:
                sub_counts = gene_counts[list(subset)]
                pan_genes = sub_counts[(sub_counts > 0).any(axis=1)].shape[0]
                core_genes = sub_counts[(sub_counts > 0).all(axis=1)].shape[0]
                pan_core_results.append({"subset": subset, "num_species": i, "pan": pan_genes, "core": core_genes})
        return pd.DataFrame(pan_core_results)

# Additional visualization functions
def plot_pie_chart(non_zero_counts, category_colors, output_dir=None, filename=None):
    """Plot pie chart of gene family classifications."""
    classifications = non_zero_counts['category'].value_counts()
    colors = [category_colors[cat] for cat in classifications.index]
    classifications.plot.pie(autopct='%1.1f%%', startangle=90, colors=colors)
    plt.ylabel("")
    plt.title("Gene Family Classifications")
    plt.savefig(os.path.join(output_dir, filename), bbox_inches='tight')
    # plt.show()
    plt.close()

def plot_core_pan_trends(core_pan_data, boxplot=False, colors=None, output_dir=None, filename=None):
    """Plot core and pan gene trends."""
    if boxplot:
        sns.boxplot(data=core_pan_data, x="num_species", y="core", color=colors["core"])
        sns.boxplot(data=core_pan_data, x="num_species", y="pan", color=colors["pan"])
    else:
        sns.scatterplot(data=core_pan_data, x="num_species", y="core", color=colors["core"])
        sns.scatterplot(data=core_pan_data, x="num_species", y="pan", color=colors["pan"])
    plt.savefig(os.path.join(output_dir, filename), bbox_inches='tight')
    # plt.show()
    plt.close()

def plot_histplot(non_zero_counts, category_colors, output_dir=None, filename=None):

    genome_number = non_zero_counts.iloc[:, 0:-1].sum(axis=1).reset_index().rename(columns={0: 'genome_number'})
    classification = non_zero_counts['category'].reset_index().rename(columns={'category': 'classification'})
    classification_count = pd.merge(genome_number, classification, on='Orthogroup', how='left').iloc[:, 1:].value_counts().reset_index(name='count')
    sorted_df = classification_count.sort_values(by='genome_number', ascending=True)

    plt.figure(figsize=(10, 6))
    sns.barplot(x='genome_number',
                y='count', hue='classification',
                data=sorted_df,
                palette=category_colors)
    plt.xlabel("Genome Number")
    plt.ylabel("Family Number")
    plt.legend(title='Proportion of pangenes')
    plt.savefig(os.path.join(output_dir, filename))
    plt.close()

def plot_heatmap(non_zero_counts, category_colors, sort=None, output_dir=None, filename=None):
    category_order = ['core', 'softcore', 'dispensable', 'private']
    other_columns = non_zero_counts.iloc[:, 0:-1].columns

    # 将分类列转为有序的分类数据类型
    non_zero_counts['category'] = pd.Categorical(non_zero_counts['category'],
                                                 categories=category_order,
                                                 ordered=True)

    if sort:
        non_zero_counts_sorted = non_zero_counts.sort_values(
            by=['category'] + list(other_columns),
            ascending=[True] + [False] * len(other_columns)
        )
    else:
        non_zero_counts_sorted = non_zero_counts.sort_values(
            by=['category']
        )

    # 将类别列中的每个值替换为颜色
    non_zero_counts_sorted['category'] = non_zero_counts_sorted['category'].apply(lambda x: category_colors[x])

    # 转换为 numpy 数组用于绘制热图
    non_zero_counts_sorted_np = non_zero_counts_sorted.iloc[:, 0:-1].T.to_numpy()
    row_labels = non_zero_counts.iloc[:, 0:-1].columns

    # 创建自定义的颜色映射
    custom_cmap = sns.color_palette(["grey", "orange"], as_cmap=True)

    # 创建热图
    fig, ax = plt.subplots(figsize=(15, 4))
    sns.heatmap(non_zero_counts_sorted_np,
                annot=False,
                cmap=custom_cmap,
                yticklabels=row_labels,
                xticklabels=False,
                ax=ax,
                cbar=True)

    # 获取热图的位置和宽度
    heatmap_pos = ax.get_position()  # 返回一个Bbox对象
    heatmap_left = heatmap_pos.xmin
    heatmap_width = heatmap_pos.width

    # 添加边际颜色条，确保宽度与热图一致
    category_numeric = non_zero_counts_sorted['category'].map(lambda x: list(category_colors.values()).index(x)).values

    # 创建边际颜色条轴，调整为热图宽度
    ax_marginal = fig.add_axes([heatmap_left, heatmap_pos.ymax + 0.02, heatmap_width, 0.05])  # 根据热图宽度调整
    marginal_cmap = ListedColormap(list(category_colors.values()))  # 使用已定义的颜色映射
    ax_marginal.imshow([category_numeric], aspect="auto", cmap=marginal_cmap)  # 用数值绘制颜色条
    ax_marginal.axis('off')

    # 手动添加图例
    legend_elements = [Patch(facecolor=color, label=label) for label, color in category_colors.items()]
    ax.legend(handles=legend_elements, loc='lower right', bbox_to_anchor=(1.05, 0), ncol=1)

    plt.title("Gene Presence/Absence Heatmap with Marginal Colors")
    plt.savefig(os.path.join(output_dir, filename), bbox_inches='tight')
    # plt.show()
    plt.close()

# Main script
def setup_parser(parser):
    parse_orthofinder_parser = parser.add_parser('parse_orthofinder', help='Analyze OrthoFinder results.')
    parse_orthofinder_parser.add_argument("-d", "--dir", required=True, help="Path to OrthoFinder result directory.")
    parse_orthofinder_parser.add_argument("-o", "--output_dir", default='.', help="Path to result directory.")
    parse_orthofinder_parser.add_argument("-p", "--pattern", required=True, choices=["pangene", "copy", "align"], help="Analysis type.")
    parse_orthofinder_parser.add_argument("--threshold_pct", type=float, default=0.9, help="pangene option: Percentage threshold for softcore genes. [0.9]")
    parse_orthofinder_parser.add_argument("--threshold_count", type=int, help="pangene option: Count threshold for softcore genes. []")
    parse_orthofinder_parser.add_argument("--max_combinations", type=int, default=1000, help="pangene option: Maximum combinations for core-pan analysis.")
    parse_orthofinder_parser.add_argument("--heatmap_sort", action="store_true", help="pangene option: Sort the heatmap. [False]")
    parse_orthofinder_parser.add_argument("--boxplot", action="store_true", help="pangene option: Plot core-pan trends as boxplot. [False]")
    parse_orthofinder_parser.add_argument("--colors", help="pangene option: Colors of gene types.")

    return parse_orthofinder_parser

def run(args):

    if args.output_dir and not os.path.exists(args.output_dir):
        os.makedirs(args.output)  # 创建输出目录

    results = OrthoFinderResults(args.dir)
    if args.pattern == "pangene":
        category_colors = {
        "core": "#d62728",
        "softcore": "#ff7f0e",
        "dispensable": "#1f77b4",
        "private": "#2ca02c",
        }
        non_zero_counts = results.classify_gene_families(threshold_pct = args.threshold_pct)
        core_pan_data = results.core_pan_analysis(1000)
        core_pan_data.to_csv(os.path.join(args.output_dir, "core_pan_analysis.csv"), index=False)
        plot_pie_chart(non_zero_counts, category_colors, output_dir=args.output_dir, filename='gene_classifications_pie.pdf')
        plot_histplot(non_zero_counts, category_colors, output_dir=args.output_dir, filename='gene_classifications_histplot.pdf')
        plot_heatmap(non_zero_counts, category_colors, args.heatmap_sort, output_dir=args.output_dir, filename='gene_classifications_heatmap.pdf')
        plot_core_pan_trends(core_pan_data, boxplot='boxplot', colors={"core": "blue", "pan": "green"}, output_dir=args.output_dir, filename='core_pan_analysis.pdf')
        df_extracted_genes = results.extract_genes(non_zero_counts)
        df_extracted_genes.to_csv(os.path.join(args.output_dir, "extracted_genes.csv"), index=False)

    elif args.pattern == "copy":
        print("Copy number analysis is not yet implemented.")
    
    else:
        print("Analysis is not yet implemented.")