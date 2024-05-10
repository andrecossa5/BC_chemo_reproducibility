"""
Plotting functions.
"""

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from circlify import circlify, Circle
from plotting_utils._plotting_base import *


## 

        
def packed_circle_plot(
    df, covariate=None, ax=None, color='b', cmap=None, alpha=.5, linewidth=1.2,
    t_cov=.01, annotate=False, fontsize=6
    ):
    """
    Circle plot. Packed.
    """
    df = df.sort_values(covariate, ascending=False)
    circles = circlify(
        df[covariate].to_list(),
        show_enclosure=True, 
        target_enclosure=Circle(x=0, y=0, r=1)
    )
    
    lim = max(
        max(
            abs(c.x) + c.r,
            abs(c.y) + c.r,
        )
        for c in circles
    )
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    
    if isinstance(color, str) and not color in df.columns:
        colors = { k : color for k in df.index }
    elif isinstance(color, str) and color in df.columns:
        colors = create_palette(df, covariate, cmap)
    else:
        assert isinstance(color, dict)
        colors = color
        print('Try to use custom colors...')

    for name, circle in zip(df.index[::-1], circles): # Don't know why, but it reverses...
        x, y, r = circle
        ax.add_patch(
            plt.Circle((x, y), r*0.95, alpha=alpha, linewidth=linewidth, 
                fill=True, edgecolor=colors[name], facecolor=colors[name])
        )

        if annotate:
            cov = df.loc[name, covariate]
            if cov > t_cov:
                n = name if len(name)<=5 else name[:5]
                ax.annotate(
                    f'{n}: {df.loc[name, covariate]:.2f}', 
                    (x,y), 
                    va='center', ha='center', fontsize=fontsize
                )

    ax.axis('off')
    
    return ax


##


def consensus_volcano(results=None, L=None, df=None, pattern=None, t=0.5, xlim=(-1,1),
                    figsize=(7,7), return_df=True, n=None):
    """
    Create a consensus volcano plot.
    """
    
    # Df
    if L is None and df is None:
        L = []
        for k in results.results:
            k_list = k.split('|')
            if bool(re.search(pattern, k_list[0])) and k_list[2] == 'wilcoxon':
                df_ = results.results[k]['df'].query('comparison == "g0_vs_g1"')
                df_.loc[:, 'comparison'] = '_'.join([k_list[0], 'vs_rest'])
                L.append(df_)
        n = len(L)
                
    elif df is None:
        
        df = (
            pd.concat(L, axis=0).reset_index()
            .rename(columns={'index':'gene'})
            .groupby('gene')['effect_size', 'evidence']
            .mean()
            .sort_values(by='effect_size', ascending=False)   
            .assign(to_annotate_positive=lambda x: (x['effect_size']>=t) & (x['evidence']<=0.1))
            .assign(to_annotate_negative=lambda x: (x['effect_size']<=-t) & (x['evidence']<=0.1))
        )
        df['log_evidence'] = -np.log10(df['evidence']+0.00000001)
        n = len(L)
        
    else:
        pass

    # Fig
    fig, ax = plt.subplots(figsize=figsize)
    scatter(df.loc[lambda x: ~x['to_annotate_positive'] & ~x['to_annotate_positive']], 
            'effect_size', 'log_evidence', 
            c='darkgrey', s=1, ax=ax)
    scatter(df.loc[lambda x: x['to_annotate_positive']], 'effect_size', 'log_evidence', 
            c='r', s=10, ax=ax)
    scatter(df.loc[lambda x: x['to_annotate_negative']], 'effect_size', 'log_evidence', 
            c='b', s=10, ax=ax)

    ax.vlines(0, df['log_evidence'].min(), df['log_evidence'].max(), 
            linestyles='dashed', colors='k')
    ax.vlines(t, df['log_evidence'].min(), df['log_evidence'].max(), 
            linestyles='dashed', colors='r')
    ax.vlines(-t, df['log_evidence'].min(), df['log_evidence'].max(),
            linestyles='dashed', colors='b')
    ax.hlines(-np.log10(0.1), -1, 1, linestyles='dashed', colors='k')
    ax.set(xlim=xlim)
    format_ax(ax, title=f'Consensus {pattern}' if pattern is not None else f'Consensus', 
            xlabel=f'Mean (n={n}) log2FC', ylabel=f'-log10(Mean (n={n}) FDR)')

    ta.allocate_text(
        fig, ax, 
        df.loc[lambda x: x['to_annotate_positive'] | x['to_annotate_negative']]['effect_size'],
        df.loc[lambda x: x['to_annotate_positive'] | x['to_annotate_negative']]['log_evidence'],
        df.loc[lambda x: x['to_annotate_positive'] | x['to_annotate_negative']].index,
        x_scatter=df['effect_size'], y_scatter=df['log_evidence'], 
        linecolor='black', textsize=8, 
        max_distance=0.5, linewidth=0.5, nbr_candidates=100
    )
    ax.spines[['top', 'right']].set_visible(False)
    
    if return_df:
        return fig, df
        
    else:
        return fig