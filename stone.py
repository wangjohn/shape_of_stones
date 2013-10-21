import random
import math
from shapely import geometry
from shapely import affinity

class VertexChipper:
    def __init__(self, stone):
        self.stone = stone

    def chip(self, iterations):
        for i in xrange(iterations):
            self.chip_vertex(self.stone)

    def chip_vertex(self, stone, index = None):
        raise "VertexChipper is an abstract class. Please implement chip_vertex"

    def random_vertex(self, polygon, index = None):
        if not index:
            index = random.randrange(len(polygon.exterior.coords)-1)

        return polygon.exterior.coords[index]

    def get_random_number(self, min_number, max_number):
        number = self.generate_random_number()
        while number < min_number or number > max_number:
            number = self.generate_random_number()

        return number

    def generate_random_number(self):
        raise "Subclasses of VertexChipper must implement this method."


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
    def __init__(self, stone, fraction_mean = 0.5, fraction_deviation = 0.15):
        self.stone = stone
        self.fraction_mean = fraction_mean
        self.fraction_deviation = fraction_deviation

    def chip_vertex(self, stone, index = None):
        vertex = self.random_vertex(stone.polygon, index)
        left_neighbor = stone.polygon.exterior.coords[index-1]
        right_neighbor = stone.polygon.exterior.coords[index+1]
        left_line = geometry.LineString([left_neighbor, vertex])
        right_line = geometry.LineString([right_neighbor, vertex])

        new_left_point = self.get_random_point_on_line(left_line)
        new_right_point = self.get_random_point_on_line(right_line)

        new_exterior = list(stone.polygon.exterior.coords)
        new_exterior[index] = new_right_point
        new_exterior.insert(index, new_left_point)

        stone.polygon = geometry.Polygon(new_exterior)

    def get_random_point_on_line(self, line):
        return line.interpolate(self.get_random_number(0.0, 1.0), True).coords[0]

    def generate_random_number(self):
        return random.normalvariate(self.fraction_mean, self.fraction_deviation)

class AngleVertexChipper(VertexChipper):
    def chip_vertex(self, stone, index = None):
        vertex = self.random_vertex(stone.polygon, index)
        centroid = stone.polygon.centroid()
        vertex_line = geometry.LineString([centroid, vertex])

        angle, is_clockwise = self.get_random_angle()
        new_line = self.create_rotated_line(base_line, angle)
        polygon_intersection = self.get_intersection(stone, line, centroid)

        stone.polygon = self.remove_vertices_between(stone.polygon, polygon_intersection, vertex, is_clockwise)

    # Removes the vertices between +point+ and +vertex+.
    def remove_vertices_between(self, polygon, point, vertex, is_clockwise=True):
        coordinates = polygon.exterior.coords
        if not is_clockwise:
            coordinates = coordinates.reverse()

        vertex_index = None
        point_index = None
        for i in xrange(len(coordinates)-1):
            current_vertex = polygon.exterior.coords[i]
            next_vertex = polygon.exterior.coords[i+1]
            if current_vertex.almost_equals(vertex):
                vertex_index = i

            line = geometry.LineString([current_vertex, next_vertex])
            if line.intersects(point):
                point_index = i
                break

        new_coordinates = coordinates[:vertex_index] + [point.coords[0]] + coordinates[point_index+1:]
        if not is_clockwise:
            new_coordinates = new_coordinates.reverse()

        return geometry.Polygon(new_coordinates)

    def get_intersection(self, stone, line, center):
        while not line.intersects(stone.polygon):
            line = affinity.scale(line, xfact=2.0, yfact=2.0, origin=center)

        return line.intersection(stone.polygon).coords[0]

    # Creates a new line which is a specified angle from the base_line.
    # The base_line should be specified so that the first point in the line is
    # the point from which the new line should intersect with the old line.
    def create_rotated_line(self, base_line, angle):
        center = base_line[0]
        return affinity.rotate(base_line, angle, center, use_radians=True)

    # Returns and angle and whether or not you are moving in the forward direction
    # in the clockwise direction
    def get_random_angle(self):
        angle = self.get_random_number(0.0, 2*math.pi)
        if angle < math.pi:
            return (angle, False)
        else:
            return (angle, True)

    def generate_random_number(self):
        return random.normalvariate(self.angle_mean, self.angle_deviation)

class Stone:
    def __init__(self, polygon):
        self.polygon = polygon
