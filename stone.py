import math

class Stone:
    def __init__(self, vertices):
        self.vertices = vertices

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

        for i in xrange(len(self.vertices)):
            j = (i+1) & (len(self.vertices) - 1)
            total_x += (self.vertices[i].x + self.vertices[j].x) * self._point_convolution(i)
            total_y += (self.vertices[i].y + self.vertices[j].y) * self._point_convolution(i)

        center_x = float(total_x) / (6.0 * signed_area)
        center_y = float(total_y) / (6.0 * signed_area)

        return Point(center_x, center_y)

    def find_signed_area(self):
        total = sum([self._point_convolution(i) for i in xrange(len(self.vertices))])
        return float(total) / 2.0

    # Returns the following:
    #
    #  x_i * y_{i+1} + x_{i+1} * y_i
    #
    def _point_convolution(self, i):
        j = (i+1) % (len(self.vertices) - 1)
        return self.vertices[i].x * self.vertices[j].y - self.vertices[j].x * self.vertices[i].y)


class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def distance(self, other):
        return math.sqrt((self.x - other.x)**2.0 + (self.y - other.y)**2.0)

    def distance_to_origin(self):
        return math.sqrt(self.x**2.0 + self.y**2.0)
