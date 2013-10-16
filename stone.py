import random
from shapely import geometry

class Stone:
    def __init__(self, polygon, fraction_mean = 0.5, fraction_deviation = 0.15):
        self.polygon = polygon
        self.fraction_mean = fraction_mean
        self.fraction_deviation = fraction_deviation

    def chip_vertex(self, index = None):
        if not index:
            index = random.randrange(len(self.polygon.exterior.coords)-1)

        vertex = self.get_vertex(index)
        left_neighbor = self.get_vertex(index-1)
        right_neighbor = self.get_vertex(index+1)
        left_line = geometry.LineString([left_neighbor, vertex])
        right_line = geometry.LineString([right_neighbor, vertex])

        new_left_point = self.get_random_point_on_line(left_line)
        new_right_point = self.get_random_point_on_line(right_line)

        new_exterior = list(self.polygon.exterior.coords)
        new_exterior[index] = new_right_point
        new_exterior.insert(index, new_left_point)

        self.polygon = geometry.Polygon(new_exterior)

    def get_vertex(self, index):
        return self.polygon.exterior.coords[index]

    def get_random_point_on_line(self, line):
        return line.interpolate(self.get_random_fraction(), True).coords[0]

    def get_random_fraction(self):
        fraction = self._generate_potential_fraction()
        while fraction < 0 or fraction > 1:
            fraction = self._generate_potential_fraction()

        return fraction

    def _generate_potential_fraction(self):
        return random.normalvariate(self.fraction_mean, self.fraction_deviation)

if __name__ == '__main__':
    from matplotlib import pyplot
    from descartes.patch import PolygonPatch

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
