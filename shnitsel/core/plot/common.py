import matplotlib as mpl


def outlabel(ax, letter):
    fixedtrans = mpl.transforms.ScaledTranslation(
        -20 / 72, +7 / 72, ax.figure.dpi_scale_trans
    )
    transform = ax.transAxes + fixedtrans
    return ax.text(
        0.0,
        1.0,
        letter,
        transform=transform,
        va='bottom',
        fontweight='bold',
        bbox=dict(facecolor='0.9', edgecolor='none', pad=3.0),
    )


def inlabel(ax, letter):
    return ax.annotate(
        letter,
        xy=(1, 1),
        xycoords='axes fraction',
        xytext=(-1, -0.5),
        textcoords='offset fontsize',
        va='top',
        fontweight='bold',
        bbox=dict(facecolor='0.9', edgecolor='none', pad=3.0),
    )
