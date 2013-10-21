import random
from shapely import geometry

class VertexChipper:
    def __init__(self, options = {}):
        self.options = options

    def chip_vertex(self, stone, index = None):
        raise "VertexChipper is an abstract class. Please implement chip_vertex"

    def random_vertex(self, polygon, index = None):
        if not index:
            index = random.randrange(len(polygon.exterior.coords)-1)

        return polygon.exterior.coords[index]

# This vertex chipper selects a random vertex to chip, then gets the lines that
# intersect vertex. It picks one point from each intersecting line and connects
# those two points to chip the stone.
#
# The points are chosen with a normal distribution centered about the
# +fraction_mean+ of the line segment. The standard deviation is given by the 
# +fraction_deviation+.
#
# A +fraction_mean+ of 0.5 means that the normal distribution is centered around
# the midpoint of the line segment. A +fraction_mean+ of 0.1 means that the normal
# distribution is centered around 0.1 times the distance of segment.
class RandomPointVertexChipper(VertexChipper):

    def chip_vertex(self, stone, index = None):
        vertex = self.random_vertex(stone.polygon, index)


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
