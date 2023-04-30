"""
Plotting functions.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from circlify import _bubbles, circlify, Circle


#

        
def packed_circle_plot(df, covariate=None, ax=None, color='b', annotate=False):
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
    
    for name, circle in zip(df.index[::-1], circles): # Don't know why, but it reverses...
        x, y, r = circle
        ax.add_patch(
            plt.Circle((x, y), r*0.95, alpha=0.5, linewidth=1.2, 
                fill=True, edgecolor=color, facecolor=color)
        )
        
        if annotate:
            cov = df.loc[name, covariate]
            if cov > 0.01:
                ax.annotate(
                    f'{name}: {df.loc[name, covariate]:.2f}', 
                    (x,y), 
                    va='center', ha='center', fontsize=5
                )

    ax.axis('off')
    
    return ax