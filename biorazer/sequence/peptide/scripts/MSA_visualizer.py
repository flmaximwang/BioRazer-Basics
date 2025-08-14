import numpy as np
from Bio import SeqIO
from matplotlib import pyplot as plt

def read_a3m(a3m_file_path):
    try:
        records = list(SeqIO.parse(a3m_file_path, "fasta"))
    except Exception as e:
        print(f"文件读取错误: {str(e)}")
        return
    
    if not records:
        print("错误: 未找到有效序列")
        return
    
    alignment = []
    for record in records:
        # 将非大写字母 (核心比对区域) 转换为 -
        cleaned_seq = "".join([c for c in str(record.seq) if not c.islower()])
        alignment.append(cleaned_seq)
    
    lengths = set([len(seq) for seq in alignment])
    if len(lengths) > 1:
        print("错误: 序列长度不一致")
        return
    
    # 将 alignment 转换为二维 numpy 数组
    return alignment

def plot_msa(msa: list[str], seq_len_list: list[int], sort_lines=True, dpi=100):
    """
    可视化多序列比对（MSA）覆盖情况（第二版）
    
    参数：
    feature_dict: 包含MSA特征的字典
    sort_lines: 是否按序列相似度排序（默认True）
    dpi: 图像分辨率
    """
    # 获取查询序列
    seq = msa.pop(0)
    seq = np.array(list(seq))
    msa = np.array([list(seq) for seq in msa])
    
    Ls = seq_len_list
    
    # 计算链的累积长度
    Ln = np.cumsum([0] + Ls)
    N = msa.shape[0]
    
    # 处理MSA数据
    gap = msa != "-"  # 检测非空位（21代表空位）
    qid = msa == seq  # 与查询序列相同的位点
    
    # 计算链级别的覆盖情况
    gapid = np.stack([gap[:,Ln[i]:Ln[i+1]].max(-1) for i in range(len(Ls))],-1)
    
    # 组织绘图数据
    lines = []
    Nn = []
    for g in np.unique(gapid, axis=0):
        i = np.where((gapid == g).all(axis=-1))
        qid_ = qid[i]
        gap_ = gap[i]
        
        # 计算序列相似度
        seqid = np.stack([qid_[:,Ln[i]:Ln[i+1]].mean(-1) for i in range(len(Ls))],-1).sum(-1) / (g.sum(-1) + 1e-8)
        
        # 处理非空位数据
        non_gaps = gap_.astype(float)
        non_gaps[non_gaps == 0] = np.nan
        
        # 排序处理
        if sort_lines:
            lines_ = non_gaps[seqid.argsort()] * seqid[seqid.argsort(),None]
        else:
            lines_ = non_gaps[::-1] * seqid[::-1,None]
        
        Nn.append(len(lines_))
        lines.append(lines_)
    
    # 合并数据
    Nn = np.cumsum(np.append(0,Nn))
    lines = np.concatenate(lines,0)
    
    # 绘制图像
    plt.figure(figsize=(8,5), dpi=dpi)
    plt.title("Sequence coverage")
    plt.imshow(lines,
              interpolation='nearest', 
              aspect='auto',
              cmap="rainbow_r", 
              vmin=0, vmax=1, 
              origin='lower',
              extent=(0, lines.shape[1], 0, lines.shape[0]))
    
    # 添加链分隔线
    for i in Ln[1:-1]:
        plt.plot([i,i],[0,lines.shape[0]],color="black")
    # 添加序列分组线
    for j in Nn[1:-1]:
        plt.plot([0,lines.shape[1]],[j,j],color="black")
    
    # 绘制覆盖率曲线
    plt.plot((np.isnan(lines) == False).sum(0), color='black')
    plt.xlim(0,lines.shape[1])
    plt.ylim(0,lines.shape[0])
    plt.colorbar(label="Sequence identity to query")
    plt.xlabel("Positions")
    plt.ylabel("Sequences")
    return plt

