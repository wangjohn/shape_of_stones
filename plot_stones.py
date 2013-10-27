from matplotlib import pyplot
from descartes.patch import PolygonPatch
from stone import *
import numpy

# Random Point Chipping
#
# fraction_mean = 0.5
# fraction_std = 0.15
# vertex_chipper = RandomPointVertexChipper(fraction_mean, fraction_std)

# Random Angle Vertex Chipping

# angle_mean = 0.35
# angle_std = 0.35
# vertex_chipper = AngleVertexChipper(angle_mean, angle_std)

# Area Chipping

area = 0.1
fraction_mean = 0.5
fraction_std = 0.15
vertex_chipper = AreaVertexChipper(area, fraction_mean, fraction_std)

def plot_coords(ax, ob):
    x, y = ob.xy
    ax.plot(x, y, 'o', color='#999999', zorder=1)

def plot_many(vertex_chipper, chips = 50):
    fig = pyplot.figure(1, dpi=90)
    for i in xrange(9):
        square = geometry.box(0.0, 0.0, 1.0, 1.0)
        stone = Stone(square)
        vertex_chipper.set_stone(stone)

        vertex_chipper.chip(chips)

        ax = fig.add_subplot(int('33' + str(i)))
        polygon = vertex_chipper.stone.polygon
        patch = PolygonPatch(polygon, facecolor='#6699cc', edgecolor='#6699cc', alpha=0.5, zorder=2)
        ax.add_patch(patch)
        ax.set_xlim(*[-0.5,1.5])
        ax.set_ylim(*[-0.5,1.5])
    pyplot.show()

def plot_single(vertex_chipper, chips = 50):
    fig = pyplot.figure(1, dpi=90)
    square = geometry.box(0.0, 0.0, 1.0, 1.0)
    stone = Stone(square)
    vertex_chipper.set_stone(stone)

    vertex_chipper.chip(chips)

    polygon = vertex_chipper.stone.polygon
    patch = PolygonPatch(polygon, facecolor='#6699cc', edgecolor='#6699cc', alpha=0.5, zorder=2)
    ax = fig.add_subplot('111')
    ax.add_patch(patch)
    ax.set_xlim(*[-0.5,1.5])
    ax.set_ylim(*[-0.5,1.5])
    pyplot.show()

def plot_histogram(vertex_chipper, iterations = 50, chips = 50):
    histogram = []
    for i in xrange(iterations):
        square = geometry.box(0.0, 0.0, 1.0, 1.0)
        stone = Stone(square)
        vertex_chipper.set_stone(stone)
        vertex_chipper.chip(chips)
        polygon = vertex_chipper.stone.polygon

        distances = vertex_chipper.stone.distances_from_centroid()
        histogram.extend(distances)

    hist = numpy.histogram(histogram)
    pyplot.hist(hist[0], hist[1])

if __name__ == '__main__':
    #plot_single(vertex_chipper, 50)
    plot_histogram(vertex_chipper)
