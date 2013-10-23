from matplotlib import pyplot
from descartes.patch import PolygonPatch
from stone import *

# Random Point Chipping
#
# fraction_mean = 0.5
# fraction_std = 0.15
# vertex_chipper = RandomPointVertexChipper(fraction_mean, fraction_std)

# Random Angle Vertex Chipping

angle_mean = 0.3
angle_std = 0.15
vertex_chipper = AngleVertexChipper(angle_mean, angle_std)

if __name__ == '__main__':
    def plot_coords(ax, ob):
        x, y = ob.xy
        ax.plot(x, y, 'o', color='#999999', zorder=1)

    fig = pyplot.figure(1, dpi=90)
    for i in xrange(9):
        square = geometry.box(0.0, 0.0, 1.0, 1.0)
        stone = Stone(square)
        vertex_chipper.set_stone(stone)

        vertex_chipper.chip(50)

        ax = fig.add_subplot(int('33' + str(i)))
        polygon = vertex_chipper.stone.polygon
        plot_coords(ax, polygon.exterior)
        patch = PolygonPatch(polygon, facecolor='#6699cc', edgecolor='#6699cc', alpha=0.5, zorder=2)
        ax.add_patch(patch)
        ax.set_xlim(*[-0.5,1.5])
        ax.set_ylim(*[-0.5,1.5])
    pyplot.show()

