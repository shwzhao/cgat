import os
import pandas as pd
from collections import defaultdict
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def classify_orthogroups(gene_count_file, softcore_threshold=0.9):
    """
    将 Orthogroups 分类为 Core、Softcore、Dispensable、Private。
    gene_count_file='Orthogroups.GeneCount.tsv'
    """
    data = pd.read_csv(gene_count_file, sep="\t")
    total_species = data.shape[1] - 2  # 排除 "Orthogroup" 和 "Total"
    data["PresentSpecies"] = (data.iloc[:, 1:-1] > 0).sum(axis=1)

    core_threshold = total_species
    softcore_threshold = int(softcore_threshold * total_species)

    classifications = {
        "Core": data[data["PresentSpecies"] == core_threshold]["Orthogroup"].tolist(),
        "Softcore": data[
            (data["PresentSpecies"] < core_threshold) & 
            (data["PresentSpecies"] >= softcore_threshold)
        ]["Orthogroup"].tolist(),
        "Dispensable": data[
            (data["PresentSpecies"] < softcore_threshold) & 
            (data["PresentSpecies"] > 1)
        ]["Orthogroup"].tolist(),
        "Private": data[data["PresentSpecies"] == 1]["Orthogroup"].tolist(),
    }
    return classifications

def extract_genes_by_species(orthogroups_file, classifications, output_file):
    """
    从 Orthogroups.tsv 中提取不同分类下的基因。
    orthogroups_file = 'Orthogroups.tsv'
    将 genes_by_category 写出为 TSV 文件。
    """
    data = pd.read_csv(orthogroups_file, sep="\t", index_col=0)
    genes_by_category = defaultdict(lambda: defaultdict(list))

    for category, orthogroups in classifications.items():
        for orthogroup in orthogroups:
            if orthogroup in data.index:
                for species, genes in data.loc[orthogroup].items():
                    if not pd.isna(genes):
                        genes_by_category[category][species].extend(genes.split(", "))
    rows = []
    for category, species_dict in genes_by_category.items():
        for species, genes in species_dict.items():
            rows.append({"Category": category, "Species": species, "Genes": "; ".join(genes)})

    df = pd.DataFrame(rows)
    df.to_csv(output_file, index=False, sep="\t")
    # print(f"结果已保存到 {output_file}")

def calculate_diversity_index(gene_counts, index_type="shannon"):
    """
    计算 Shannon 或 Simpson 多样性指数。
    """
    total_genes = sum(gene_counts)
    proportions = [count / total_genes for count in gene_counts if count > 0]

    if index_type == "shannon":
        return -sum(p * np.log(p) for p in proportions)
    elif index_type == "simpson":
        return 1 - sum(p ** 2 for p in proportions)
    else:
        raise ValueError("Invalid index type. Choose 'shannon' or 'simpson'.")


def shared_and_pangenes(data, max_combos=1000):
    """
    统计共享基因家族数量和 Pan-gene 数量的趋势。
    """
    species = data.columns[1:-1]
    combo_counts = []
    for n in range(2, len(species) + 1):
        combos = list(combinations(species, n))[:max_combos]  # 限制组合数量
        for combo in combos:
            # 共享基因家族数量
            shared_count = data[(data[list(combo)] > 0).all(axis=1)].shape[0]
            # Pan-gene 基因家族并集数量
            pan_count = data[(data[list(combo)] > 0).any(axis=1)].shape[0]
            combo_counts.append((len(combo), shared_count, pan_count))
    return pd.DataFrame(combo_counts, columns=["SpeciesCount", "SharedFamilies", "PanGenes"])

def visualize_results(classifications, shared_df, diversity_indices, output_dir):
    """
    可视化分析结果，包括分类饼图、共享基因趋势图、热图。
    """
    os.makedirs(output_dir, exist_ok=True)

    # 1. 分类饼图
    counts = [len(v) for v in classifications.values()]
    labels = classifications.keys()
    plt.pie(counts, labels=labels, autopct='%1.1f%%')
    plt.title("Orthogroup Classification")
    plt.savefig(f"{output_dir}/classification_piechart.png")
    plt.clf()

    # 2. 创建共享基因家族和 Pan-gene 的趋势箱线图
    melted_df = shared_df.melt(id_vars=["SpeciesCount"], 
                               value_vars=["SharedFamilies", "PanGenes"],
                               var_name="Category", 
                               value_name="GeneFamilies")

    plt.figure(figsize=(10, 6))
    sns.boxplot(data=melted_df, x="SpeciesCount", y="GeneFamilies", hue="Category",
                palette={"SharedFamilies": "blue", "PanGenes": "orange"})

    # 绘制拟合曲线
    for category, color in zip(["SharedFamilies", "PanGenes"], ["blue", "orange"]):
        subset = shared_df[["SpeciesCount", category]].rename(columns={category: "GeneFamilies"})
        subset = subset.sort_values("SpeciesCount")  # 按 SpeciesCount 排序

        # 拟合曲线
        x = subset["SpeciesCount"]-2
        y = subset["GeneFamilies"]
        z = np.polyfit(x, y, 2)  # 二次多项式拟合
        p = np.poly1d(z)

        # 绘制拟合曲线
        plt.plot(x, p(x), color=color, label=f"{category} Trend (Fitted)")

    plt.title("Shared and Pan-gene Families by Species Count")
    plt.xlabel("Number of Species")
    plt.ylabel("Number of Gene Families")
    plt.legend(title="Category")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/shared_and_pangenes_trend.png")
    plt.clf()

    # 3. 基因多样性热图
    diversity_df = pd.DataFrame(diversity_indices).T
    sns.heatmap(diversity_df, annot=True, cmap="coolwarm", cbar_kws={'label': 'Diversity Index'})
    plt.title("Species Diversity Indices")
    plt.savefig(f"{output_dir}/diversity_heatmap.png")
    plt.clf()


def pangene(input_dir, output_dir, softcore_threshold=0.9, max_combos=1000):
    '''
    主函数：分析多个 OrthoFinder 结果目录。
    Orthogroups/Orthogroups.GeneCount.tsv
    '''
    os.makedirs(output_dir, exist_ok=True)
    diversity_indices = {}

    # for subdir in os.listdir(input_dir):
    #     result_dir = os.path.join(input_dir, subdir)
    gene_count_file = os.path.join(input_dir, "Orthogroups.GeneCount.tsv")
    orthogroups_file = os.path.join(input_dir, "Orthogroups.tsv")

    if not os.path.exists(gene_count_file) or not os.path.exists(orthogroups_file):
        # print(f"跳过目录 {result_dir}，缺少必要文件。")
        # break
        print("file does exist.")
        import sys
        sys.exit()

    # 分类基因家族
    classifications = classify_orthogroups(gene_count_file, softcore_threshold)

    # 提取基因名称
    extract_genes_by_species(orthogroups_file, classifications, os.path.join(output_dir, "genes_by_category.tsv"))

    # 计算多样性指数
    gene_counts = pd.read_csv(gene_count_file, sep="\t").iloc[:, 1:-1]
    for species in gene_counts.columns:
        diversity_indices[species] = {
            "Shannon": calculate_diversity_index(gene_counts[species].tolist(), "shannon"),
            "Simpson": calculate_diversity_index(gene_counts[species].tolist(), "simpson")
        }

    # 统计共享基因家族和 Pan-gene 的趋势
    data = pd.read_csv(gene_count_file, sep="\t")
    shared_df = shared_and_pangenes(data, max_combos)

    # 可视化
    visualize_results(classifications, shared_df, diversity_indices, output_dir)









def parse_align():
    '''
    根据以下文件信息，获取比对结果，
    # WorkingDirectory/SpeciesIDs.txt
    # WorkingDirectory/Species*.fa
    # WorkingDirectory/Blast3_2.txt
    '''
    pass

# paralog

def retree():
    '''
    根据单拷贝OG，重新选择策略进行重新建树
    # fasttree, iq-tree, raxml, astral (https://github.com/smirarab/ASTRAL/)

    需要学习掌握直系同源和旁系同源概念。
    '''
    pass










def setup_parser(parser):
    parse_orthofinder_parser = parser.add_parser('parse_orthofinder', help='Parse orthofinder\'s results.')
    parse_orthofinder_parser.add_argument('-p', '--pattern', required=True, help='{align, retree, pangene}. pangene: Orthogroup-based pangenome analysis.')
    parse_orthofinder_parser.add_argument('-f', '--orthofinder_result_file', required=True, help='Path to orthofinder result file')
    parse_orthofinder_parser.add_argument("-o", "--output_file", required=True, help="Directory for output results.")
    parse_orthofinder_parser.add_argument("--softcore", type=float, default=0.9, help="Proportion for softcore (default: 0.9).")
    parse_orthofinder_parser.add_argument("--max-combos", type=int, default=1000, help="Max number of species combinations to analyze.")

    return parse_orthofinder_parser


def run(args):
    if args.pattern == "pangene":
        input_dir = args.orthofinder_result_file + "/Orthogroups"  # 输入目录
        # output_dir = "path_to_output"
        pangene(input_dir, args.output_file, softcore_threshold=0.9, max_combos=1000)
    else:
        pass
