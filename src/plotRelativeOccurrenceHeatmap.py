import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.colors as mcolors
import matplotlib as mpl
import textdistance as td
import sklearn.cluster
from datetime import date as dt
import sys
import json
from webcolors import rgb_to_hex
import math
mpl.rcParams['pdf.fonttype'] = 42


def get_name(f):
    return f.split('/')[-1].split('_')[0]

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    '''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    '''
    return [v/256 for v in value]

def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp


def get_unique_string_tokens_w_occurrence(test_str):
    """
    Tokenizes k-mer into all possible substrings and assigns
    order to tokens. For example ‘AGGU’ is tokenized into:
    'A1', 'AG1', 'AGG1', 'AGGU1', 'G1', 'GG1', 'GGU1', 'G2', 'GU1', 'U1'.
    It generates tokens in the RNA alphabet (ACGU), as well as
    purine - pyrimidine tokens.
    """
    tokens = [test_str[i: j] for i in range(len(test_str)) for j in range(i + 1, len(test_str) + 1)]
    RY_signature = test_str.replace('U', 'Y').replace('A', 'R').replace('G', 'R').replace('C', 'Y')
    RY_tokens = [RY_signature[i: j] for i in range(len(test_str)) for j in range(i + 1, len(test_str) + 1)]
    ugca_tokens = []
    ry_tokens = []
    unique_tokens = sorted(list(set(tokens)))
    token_counts = [tokens.count(i) for i in unique_tokens]
    ry_unique = sorted(list(set(RY_tokens)))
    ry_token_counts = [RY_tokens.count(i) for i in ry_unique]
    for tup in list(zip(unique_tokens, token_counts)):
        for i in range(tup[1]):
            ugca_tokens.append(f'{tup[0]}{i+1}')
    for tup in list(zip(ry_unique, ry_token_counts)):
        for i in range(tup[1]):
            ry_tokens.append(f'{tup[0]}{i+1}')
    return ugca_tokens, ry_tokens


def get_distances(words):
    """
    Computes pairwise jaccard distance between tokenized k-mers.
    """
    words = np.asarray(words) #So that indexing with a list will work
    similarity_ugca = np.array([[td.jaccard.similarity(get_unique_string_tokens_w_occurrence(w1)[0], get_unique_string_tokens_w_occurrence(w2)[0]) for w1 in words] for w2 in words])
    similarity_ry = np.array([[td.jaccard.similarity(get_unique_string_tokens_w_occurrence(w1)[1], get_unique_string_tokens_w_occurrence(w2)[1]) for w1 in words] for w2 in words])
    return similarity_ugca, similarity_ry


def combine_distances(dist1, dist2, x2):
    """
    Combines distance matrices.
    """
    sum_test = np.add(dist1, dist2 * x2)
    div_test = sum_test / (1 + x2)
    return div_test


def get_clustering(similarity_matrix, x2, d):
    """
    Clusters k-mers with affinity propagation, based on precomputed
    similarity matrix.
    """
    affprop = sklearn.cluster.AffinityPropagation(affinity="precomputed", damping=d, max_iter=1000, convergence_iter=200)
    affprop.fit(similarity_matrix)
    return affprop.labels_, len(np.unique(affprop.labels_))

def get_clustered_df(motifs, d, a):
    """
    Clusters k-mers by sequence and returns them as dataframe.
    """
    df_cl = pd.DataFrame(columns=['motif', 'cluster'])
    df_cl['motif'] = motifs
    similarity_ugca, similarity_ry = get_distances(motifs)
    combined = combine_distances(similarity_ugca, similarity_ry, a)
    labels, no_clusters = get_clustering(combined, a, d)
    df_cl['cluster'] = labels
    df_cl.sort_values(by='cluster', ascending=True, inplace=True)
    return df_cl

def smooth_df(df_smooth, smoot):
    """
    Smooths values in rtxn table.
    """
    new_id = [int(i) for i in df_smooth.index]
    df_smooth.set_index([new_id], drop=True, inplace=True)
    df_smooth = df_smooth.rolling(smoot, center=True, win_type="triang").mean()
    # slicing drops edge values that get NaN due to rolling mean
    df_smooth = df_smooth.iloc[int(smoot / 2) : -(int(smoot / 2) + 1), :]
    return df_smooth

def get_hlines(df_cl, cl_order):
    """
    Gets location of horizontal lines that separate k-mer clusters in heatmap.
    """
    hlines = [0]
    for cl in cl_order:
        count = df_cl['cluster'].values.tolist().count(cl)
        hlines.append(hlines[-1]+count)
    return hlines

def get_ordered_motifs(df_cl, cl_order):
    """
    Orders k-mers for heatmap visualisation.
    Within clusters, k-mers are ordered by decreasing PEKA-score.
    """
    ordered_motifs = []
    for cl in cl_order:
        df_t = df_cl.loc[df_cl['cluster'] == cl].sort_values(by='PEKA-score', ascending=False)
        ordered_motifs.extend(df_t['motif'].values.tolist())
    return ordered_motifs

def create_labels(df_cl, ordered_motifs):
    """
    Pads k-mers with underscores relative to maxpeak position
    to generate heatmap labels.
    """
    df_t = df_cl.copy()
    labels = []
    kmer_length = len(ordered_motifs[0])
    for m in ordered_motifs:
        maxp = df_t.loc[df_t['motif'] == m, 'maxp'].values[0]
        # If maxpeak is within 3 nt away from the crosslink site, pad the label with underscores
        if abs(maxp) <= 3:
            pad = (abs(maxp) * '_')
            if maxp < 0:
                #  If peak position  is upstream from crosslink, pad on the right of the motif
                l = (((kmer_length + 2*3)-len('___' + pad + m))* '_') + m + pad + '___'
            elif maxp > 0:
                #  If peak position  is dowstream from crosslink, pad on the left of the motif
                l = '___' + pad + m + (((kmer_length + 2*3)-len('___' + pad + m))* '_')
            elif maxp == 0:
                l = '___' + m + '___'
        else:
            # If maxpeak is more than 3 nt away from the crosslink site, pad the label with '...'
            pad = '...'
            if maxp < 0:
                l = m + pad
            elif maxp > 0:
                l = pad + m
        labels.append(l)
    return labels

def get_rgba(v, norm, c):
    cmap1 = mpl.cm.get_cmap(c)
    return cmap1(norm(v),bytes=True)[:-1]

def get_representative_colors(c):
    cmap = mpl.cm.get_cmap(c)
    clr = [rgb_to_hex(cmap(v, bytes=True)[:-1]) for v in list(np.linspace(0.0, 1.0, num=11))]
    return clr

################################################################################
# To save json files

def get_matrix(df, norm, c, face_col='#f6ff6d'):
    matrix = df.values
    new_matrix = []
    for el in matrix:
        new_el = []
        for v in el:
            cell = {}
            if math.isnan(v):
                cell['value'] = None
                cell['color'] = face_col
            else:
                cell['value'] = v
                cell['color'] = rgb_to_hex(get_rgba(v, norm, c))
            new_el.append(cell)
        new_matrix.append(new_el)
    return new_matrix

def get_full_dict(df, norm, c_list, hlines, vlines, xlim, ylim, 
                  peka_df, peka_cmap, peka_norm):
    dict_out = {'rbp_heatmap': {}, 'PEKA_score_heatmap' : {}}
    # For RBP heatmap save this
    cmap = get_continuous_cmap(c_list)
    dict_out['rbp_heatmap']['colors'] = c_list
    dict_out['rbp_heatmap']['columns'] = df.columns.tolist()
    dict_out['rbp_heatmap']['rows'] = df.index.tolist()
    dict_out['rbp_heatmap']['matrix'] = get_matrix(df, norm, cmap)
    dict_out['rbp_heatmap']['hlines'] = {'line_positions': hlines, 'line_start_end': xlim, 'linestyles': 'solid', 'colors': 'white'}
    dict_out['rbp_heatmap']['vlines'] = {'line_positions': vlines, 'line_start_end': ylim, 'linestyles': 'dashed', 'colors': 'white'}
    dict_out['rbp_heatmap']['yticks_font'] = 'FreeMono'
    # For PEKA-score supplementary heatmap save this
    dict_out['PEKA_score_heatmap']['cmap'] = get_representative_colors(peka_cmap)
    dict_out['PEKA_score_heatmap']['columns'] = peka_df.columns.tolist()
    dict_out['PEKA_score_heatmap']['rows'] = peka_df.index.tolist()
    dict_out['PEKA_score_heatmap']['matrix'] = get_matrix(peka_df, peka_norm, peka_cmap)
    return dict_out

def np_encoder(object):
    if isinstance(object, np.generic):
        return object.item()

################################################################################

mpl.rcParams["font.monospace"] = "Courier"


def get_figure(f, rtxn_f, output_path, w=25, n_top=40, s="kmer_length", SaveJson=False):
    date = dt.today().strftime("%Y%m%d")
    pal = ["#54478c","#2c699a","#048ba8","#0db39e","#16db93","#83e377","#b9e769","#efea5a","#f1c453","#f29e4c"]

    prot = get_name(f)
    print(prot)
    print('PEKA-score file:', f)
    print('File with relative occurrences:', rtxn_f)
    df_data = pd.read_csv(f, sep='\t', index_col=0)
    #Get motifs for clustering.
    top_m = df_data.sort_values(by='PEKA-score', ascending=False).head(n_top).index.to_list()

    if s=='kmer_length':
        s = len(top_m[0])
        print('smoothing window set to', s)
    else:
        print('smoothing window set to', s)
    #Cluster motifs based on sequence
    df_cl = get_clustered_df(top_m, 0.5, 0)
    #Import and smooth RTXN data for plotting
    df_rtxn = pd.read_csv(rtxn_f, sep='\t', index_col=0).T
    df_smooth = smooth_df(df_rtxn, s)
    df_cl['maxp'] = df_smooth.loc[-w:w, df_cl['motif'].values.tolist()].idxmax().values.tolist()
    #Order clusters based on max PEKA-score within cluster
    df_cl['PEKA-score'] = df_data.loc[df_cl['motif'].values.tolist(), 'PEKA-score'].values.tolist()
    # Sort clusters by max PEKA-score
    df_order = df_cl.groupby('cluster', as_index=False)['PEKA-score'].max().sort_values(by='PEKA-score', ascending=False)
    cl_order = df_order['cluster'].values.tolist()
    #Get horizontal lines for plotting and order motifs for heatmap
    hlines = get_hlines(df_cl, cl_order)
    #Order motifs within cluster based on their PEKA-score
    ordered_mots = get_ordered_motifs(df_cl, cl_order)
    labels = create_labels(df_cl, ordered_mots)
    df_heatmap = df_smooth.loc[-w:w, ordered_mots].T.copy()
    df_heatmap['idx'] = labels
    df_heatmap.set_index('idx', inplace=True)
    #Plotting
    vlines = [w, w+1]
    fig = plt.figure(figsize=(17, 20))
    gridspec_kw={'width_ratios': [20, 1, 0.5],
                  'wspace': 0.06}
    gs = fig.add_gridspec(2, 3, **gridspec_kw)
    ax = fig.add_subplot(gs[:2, 0])
    vmin = df_heatmap.min().min()
    vmax = df_heatmap.max().max()
    rbp_minmax = {'vmin':vmin, 'vmax':vmax}
    g = sns.heatmap(df_heatmap, ax=ax, cmap=get_continuous_cmap(pal), xticklabels=True, yticklabels=True, cbar=False, vmin=vmin, vmax=vmax)
    for _, spine in g.spines.items():
        spine.set_visible(True)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.hlines(hlines, *ax.get_xlim(), colors='white')
    ax.vlines(vlines, *ax.get_ylim(), colors='white', linestyles='dashed')
    ax.tick_params(left=False)
    ax.tick_params(axis='x', labelrotation=90)
    ax.set_ylabel('')
    ax.set_yticklabels(ax.get_yticklabels(), fontname='FreeMono', fontsize=16)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=14)
    ax.set_xlabel('nt position relative to tXn', fontsize=16, labelpad=1)
    ax.set_title(f'{prot} - Relative occurence around tXn for top {n_top} motifs', fontsize=18)
    cax = fig.add_subplot(gs[0, 2])
    cax.tick_params(labelsize=16)
    cb_norm_heatmap = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    fig.colorbar(mpl.cm.ScalarMappable(norm=cb_norm_heatmap, cmap=get_continuous_cmap(pal)), cax=cax, ax=ax)
    # cb_ticks = list(cax.get_yticks())
    # Get only visible ticks
    x0, x1 = cax.get_ylim()
    cb_ticks = [t for t in cax.get_yticks() if t>=x0 and t<=x1]
    # Plot PEKA-score
    ax = fig.add_subplot(gs[:2, 1])
    vmin = df_cl['PEKA-score'].min()
    vmax = df_cl['PEKA-score'].max()
    peka_minmax = {'vmin':vmin, 'vmax':vmax}
    df_peka = df_data.loc[ordered_mots, 'PEKA-score'].to_frame()
    df_peka = df_peka.rename(index=dict(zip(ordered_mots, labels)))
    g = sns.heatmap(df_peka, ax=ax, vmin=vmin, vmax=vmax, cmap='gist_yarg', xticklabels=True, yticklabels=False, cbar=False)
    for _, spine in g.spines.items():
        spine.set_visible(True)
    ax.set_xticklabels(labels=['PEKA-score'], fontsize=16, rotation=-30, ha='left')
    ax.tick_params(bottom=True)
    cax = fig.add_subplot(gs[1, 2])
    cb_norm_peka = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    fig.colorbar(mpl.cm.ScalarMappable(norm=cb_norm_peka, cmap='gist_yarg'), cax=cax, ax=ax)
    cax.tick_params(labelsize=16)
    # peka_ticks = list(cax.get_yticks())
    # Get only visible ticks
    x0, x1 = cax.get_ylim()
    peka_ticks = [t for t in cax.get_yticks() if t>=x0 and t<=x1]
    plt.close()
    fig.savefig(f'{output_path}/{prot}_heatmap_top{n_top}_clustered_rtxn_underscores.pdf', bbox_inches='tight')
    if SaveJson:
        #Save dict
        d1 = get_full_dict(df_heatmap, cb_norm_heatmap, pal, hlines, vlines, xlim, ylim,
                        df_peka, 'gist_yarg', cb_norm_peka)
        d1['rbp_heatmap']['colorbar_ticks'] = [round(v, 2) for v in cb_ticks]
        d1['rbp_heatmap']['colorbar_vmin_vmax'] = rbp_minmax
        d1['PEKA_score_heatmap']['colorbar_ticks'] = [round(v, 2) for v in peka_ticks]
        d1['PEKA_score_heatmap']['colorbar_vmin_vmax'] = peka_minmax
        with open(f'{output_path}/{prot}.json', 'w') as fp:
            json.dump(d1, fp, default=np_encoder, sort_keys=True, indent=4)

def main():
    f = sys.argv[1]
    rtxn_f = sys.argv[2]
    output_path = sys.argv[3]
    w = int(sys.argv[4])
    n_top = int(sys.argv[5])

    get_figure(f,
    rtxn_f,
    output_path,
    w=w,
    n_top=n_top
    )

if __name__=='__main__':
    main()

