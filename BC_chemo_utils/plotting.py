"""
Plotting functions.
"""

import matplotlib.pyplot as plt
import circlify


##


def packed_circle_plot(df, covariate=None, ax=None, color='b', annotate=False):
    """
    Circle plot. Packed.
    """
    circles = circlify.circlify(
        df[covariate].tolist(), 
        show_enclosure=True, 
        target_enclosure=circlify.Circle(x=0, y=0, r=1)
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

    for name, circle in zip(df.index, circles):
        x, y, r = circle
        ax.add_patch(
            plt.Circle((x, y), r*0.95, alpha=0.4, linewidth=1, 
                fill=True, edgecolor="black", facecolor=color)
        )
        if annotate:
            ax.annotate(name, (x,y), va='center', ha='center', fontsize=4)
    ax.axis('off')

    return ax