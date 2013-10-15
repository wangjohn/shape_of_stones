import math

class Polygon:
    # Vertices should be included in a clockwise order along the polygon. Vertices
    # next to each other in the list will be assumed to have edges connecting them.
    def __init__(self, vertices):
        self.vertices = vertices
        self.n = len(vertices)

    def get_vertex(self, index):
        i = (index % self.n)
        return self.vertices[i]

    # Find the centroid of the stone based on the vertices. In other words, it is
    # the center of mass assuming uniform density.
    #
    # Formula for computing the centroid is given on Wikipedia at 
    # http://en.wikipedia.org/wiki/Centroid.
    #
    def find_centroid(self):
        signed_area = self.find_signed_area()

        total_x = 0
        total_y = 0

        for i in xrange(self.n):
            total_x += (self.get_vertex(i).x + self.get_vertex(i+1).x) * self._point_convolution(i)
            total_y += (self.get_vertex(i).y + self.get_vertex(i+1).y) * self._point_convolution(i)

        center_x = float(total_x) / (6.0 * signed_area)
        center_y = float(total_y) / (6.0 * signed_area)

        return Point(center_x, center_y)

    def find_signed_area(self):
        total = sum([self._point_convolution(i) for i in xrange(self.n)])
        return float(total) / 2.0

    # Returns the following:
    #
    #  x_i * y_{i+1} + x_{i+1} * y_i
    #
    def _point_convolution(self, i):
        return self.get_vertex(i).x * self.get_vertex(i+1).y - self.get_vertex(i+1).x * self.get_vertex(i).y)


class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def distance(self, other):
        return math.sqrt((self.x - other.x)**2.0 + (self.y - other.y)**2.0)

    def distance_to_origin(self):
        return math.sqrt(self.x**2.0 + self.y**2.0)
