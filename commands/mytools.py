import sys
import numpy as np
import pandas as pd
import re

def read_data(args):
    '''
    读取基因表达矩阵，可以管道输入，也可以文件名输入，直接读取，默认有表头
    '''
    if args.input == "-":
        input = sys.stdin
    else:
        input = args.input

    # if args.header:
    #     header = 0
    # else:
    #     header = None

    data = pd.read_csv(input, sep=args.sep, index_col=0, header=0)

    return data


class GeneExpressionData:
    '''
    基因表达矩阵类
    '''
    def __init__(self, data):
        self.data = data


    def get_data(self):
        return self.data


    def get_genes(self):
        data = self.get_data()
        genes = list(data.index)
        return genes


    def get_samples(self):
        data = self.get_data()
        samples = [str(column) for column in data.columns]
        return samples


    def get_groups(self, sample_sep):
        samples = self.get_samples()
        if sample_sep == "":
            groups = [re.sub(r"\d+$", "", sample) for sample in samples]
            repeats = [re.findall(r'\d+', sample)[-1] for sample in samples]
        else:
            groups = [re.split(sample_sep, sample, 1)[0] for sample in samples]
            repeats = [re.split(sample_sep, sample, 1)[1] for sample in samples]
        return groups,repeats


    def transpose(self):
        '''
        转置
        '''
        data = self.get_data()
        index_name = data.index.name
        data_transposed = data.T.reset_index()
        data_transposed = data_transposed.rename(columns={'index': index_name})
        return data_transposed


    def longer(self, var_name="Samples", value_name="Values"):
        '''
        宽边长
        '''
        data = self.get_data()
        data_columns = data.columns
        index_name = data.index.name
        data_longered = pd.melt(data.reset_index(), id_vars=[index_name], value_vars=data_columns, var_name=var_name, value_name=value_name).sort_values(by=[index_name, var_name])
        return index_name, data_longered


    def wider(self, columns="Samples", values="Values"):
        '''
        宽边长
        '''
        data = self.get_data()
        data_coluns = data.columns
        index_name = data.index.name
        data_widered = pd.pivot_table(data.reset_index(), index=index_name, columns=columns, values=values).reset_index()
        return data_widered


    def mean(self, sample_sep):
        '''
        求均值
        '''
        samples = self.get_samples()
        groups,repeats = self.get_groups(sample_sep=sample_sep)
        group_dict = dict(zip(samples, groups))
        index_name, data_longered = self.longer()
        data_longered['Groups'] = data_longered['Samples'].map(group_dict)
        data_longered_mean = data_longered.groupby(by=[index_name, 'Groups']).mean()
        data_longered_mean_widered = pd.pivot_table(data_longered_mean, index=index_name, columns='Groups', values='Values').reset_index()
        return data_longered_mean_widered


    def round(self, digit):
        '''
        表达量的位数
        '''
        data = self.get_data()
        data_rounded = data.round(decimals=int(digit)).reset_index()
        return data_rounded


    def log(self, base):
        '''
        对表达矩阵取log, 要有识别e的能力
        '''
        # import math
        data = self.get_data()
        data_logged = np.log(self.data).reset_index()
        return data_logged


    def scaler_minmax(self, max=1, min=0):
        '''
        min-max标准化 (归一化), 数据缩放为[0, 1]之间, 然后再变换范围
        '''
        from sklearn.preprocessing import MinMaxScaler
        scaler = MinMaxScaler()
        data_scaled = scaler.fit_transform(self.data)
        data_scaled_df = pd.DataFrame(data_scaled, columns=data_scaled.columns, index=data_scaled.index)
        return data_scaled_df


    def scaler_z_score(self):
        '''
        z-score标准化
        '''
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler()
        data_standardized = scaler.transform(self.data)
        data_standardized_df = pd.DataFrame(data_standardized, columns=data_standardized.columns, index=data_standardized.index)
        return data_standardized_df


    def count2tpm(self, exon_length):
        '''
        将gff也做成一个类, 设置一些功能, 这里可以导入读取exon长度
        '''
        return None


    def filter_patten(self, patten, regex=False):
        '''
        根据命令行中输入的基因id进行筛选, -r时支持正则匹配
        '''
        return None


    def filter_geneid(self, geneid):
        '''
        根据含有基因名的文件筛选
        '''
        return None

    def filter_group(self, group=1, sample=1, exp=0.5, sep="_"):
        '''
        根据分组中表达样本的数量筛选，例如`group=1, sample=2, exp=0.5`意味着有一个组中的两个样本的表达量超过0.5, 视为这个基因表达。需要有分组的判定。
        '''
        return None


    def filter_mean(self, group=1, mean=0.5, sep="_"):
        '''
        根据分组的平均值筛选，例如`group=1, mean=0.5`表示某个基因只要有一个分组的平均值大于0.5视为表达。
        '''
        return None

    def select(self, groups, others=False, sep="_"):
        '''
        根据组名进行选取，`others=True`时，把剩下的列按文件原来的顺序依此排列到后面。
        '''
        return None


    def parse_columns_arg(self, cols):
        data = self.get_data()
        col_numbers = list(range(len(data.columns)+1))

        selected_cols = []
        ranges = cols.split(',')
        for r in ranges:
            if '-' in r:
                # start, end = map(int, r.split('-'))
                start, end = r.split('-')
                start = int(start)
                if end == "":
                    end = int(len(data.columns)+1)
                else:
                    end = int(end)
                selected_cols += range(start, end+1)
            else:
                selected_cols.append(int(r))

        selected_cols = list(map(lambda x: x - 1, selected_cols))
        other_cols = [x for x in col_numbers if x not in selected_cols]

        return data, selected_cols, other_cols


    def cut(self, cols, others=False):
        '''
        根据列号进行选取, 模仿shell中cut的参数识别和功能, 区别就是可重复选择, 可调换顺序
        '''
        data, selected_cols, other_cols = self.parse_columns_arg(cols=cols)

        df = data.reset_index()

        if others:
            df_all = pd.concat([df.iloc[:, selected_cols], df.iloc[:, other_cols]], axis=1)
        else:
            df_all = df.iloc[:, selected_cols]

        return df_all


    def bind(self, files, by):
        '''
        数据链接，有很多情况，需要仔细考虑
        '''
        return None


    def pca(self):
        '''
        主成分分析
        '''
        return None


    def cor(self, method="pearson"):
        '''
        相关性分析
        '''
        return None


    def cluster(self):
        '''
        聚类分析
        '''
        return None


    def differentexpres(self):
        '''
        差异表达分析, 引入R包
        '''
        return None


    def enricher(self, geneset):
        '''
        富集分析, 只针对非模式物种, 需要提供基因的GO注释, 背景集
        '''
        return None
