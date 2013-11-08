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

        distances = vertex_chipper.stone.distances_from_centroid(300)
        histogram.extend(distances)

    graph_histogram(histogram)

def graph_histogram(histogram):
    ax = pyplot.axes()
    hist = numpy.histogram(histogram, [i*0.05 for i in xrange(18)])
    pos = numpy.arange(len(hist[0]))
    ax.set_xticks(pos + 1.0/2)
    ax.set_xticklabels(hist[1])
    ax.set_xlabel('Scaled Distance from Centroid')
    pyplot.bar(pos, hist[0], 1.0, color='r')
    pyplot.show()

def plot_ellipse(a, b):
    histogram = []
    num_points = 20000
    max_dist = max(a, b)
    for i in xrange(num_points):
        theta = random.uniform(-math.pi / 2.0, math.pi / 2.0)
        slope = math.tan(theta)
        denominator = math.sqrt(1.0 / a + (slope**2.0) / b)
        x = 1.0 / denominator
        ysq = float(b) - (float(b) / float(a))*(x**2)
        distance = (x**2.0 + ysq)
        histogram.append(distance / max_dist)

    max_dist = max(histogram)
    min_dist = min(histogram)
    histogram = [(i - min_dist) / (max_dist - min_dist) for i in histogram]
    graph_histogram(histogram)

def plot_kite():
    kite = geometry.Polygon([(1.0, 0.0), (0.5, 1.3), (0.0, 2.5), (-0.5, 1.3), (-1.0, 0.0), (0.0, -1.0), (1.0, 0.0)])
    stone = Stone(kite)
    histogram = stone.distances_from_centroid(4000)
    print kite.centroid

    max_dist = max(histogram)
    print max_dist
    histogram = [float(i) / max_dist for i in histogram]

    graph_histogram(histogram)

if __name__ == '__main__':
    #plot_single(vertex_chipper, 50)
    plot_histogram(vertex_chipper, iterations = 100)
    #plot_many(vertex_chipper, 50)
    #plot_ellipse(1.0, 1.5)
    #plot_kite()

