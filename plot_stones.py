from matplotlib import pyplot
from descartes.patch import PolygonPatch
from stone import *

if __name__ == '__main__':
    def plot_coords(ax, ob):
        x, y = ob.xy
        ax.plot(x, y, 'o', color='#999999', zorder=1)

    fig = pyplot.figure(1, dpi=90)
    for i in xrange(9):
        square = geometry.box(0.0, 0.0, 1.0, 1.0)
        stone = Stone(square, 0.1, 0.15)

        for j in xrange(50):
            stone.chip_vertex()

        ax = fig.add_subplot(int('33' + str(i)))
        plot_coords(ax, stone.polygon.exterior)
        polygon = stone.polygon
        patch = PolygonPatch(polygon, facecolor='#6699cc', edgecolor='#6699cc', alpha=0.5, zorder=2)
        ax.add_patch(patch)
        ax.set_xlim(*[-0.5,1.5])
        ax.set_ylim(*[-0.5,1.5])
    pyplot.show()

