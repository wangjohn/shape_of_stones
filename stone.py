import random
import math
from shapely import geometry
from shapely import affinity

class VertexChipper:
    def __init__(self):
        pass

    def set_stone(self, stone):
        self.stone = stone

    def chip(self, iterations):
        for i in xrange(iterations):
            result = self.chip_vertex(self.stone)
            self.stone = result

    def chip_vertex(self, stone, index = None):
        raise "VertexChipper is an abstract class. Please implement chip_vertex"

    def random_vertex(self, polygon, index = None):
        if not index:
            index = random.randrange(len(polygon.exterior.coords)-1)

        return (polygon.exterior.coords[index], index)

    def get_intersection(self, vertex_line, line, center):
        while not line.crosses(vertex_line):
            line = affinity.scale(line, xfact=2.0, yfact=2.0, origin=center)

        object_of_intersection = line.intersection(vertex_line)
        for coord in object_of_intersection.coords:
            if not self.almost_equals(center, coord):
                return geometry.Point(coord)

    def get_random_number(self, min_number, max_number):
        number = self.generate_random_number()
        while number < min_number or number > max_number:
            number = self.generate_random_number()

        return number

    def create_line_from_center(self, angle, center):
        cx, cy = center
        line = geometry.LineString([center, (cx + 1.0, cy)])
        return self.create_rotated_line(line, angle, center)

    def generate_random_number(self):
        raise "Subclasses of VertexChipper must implement this method."

class AreaVertexChipper(VertexChipper):
    def __init__(self, area = 0.1, fraction_mean = 0.5, fraction_std = 0.15):
        self.area = area
        self.fraction_mean = fraction_mean
        self.fraction_std = fraction_std

    def chip_vertex(self, stone, index = None):
        vertex, index = self.random_vertex(stone.polygon, index)
        left_neighbor = stone.polygon.exterior.coords[index-1]
        right_neighbor = stone.polygon.exterior.coords[index+1]
        left_line = geometry.LineString([left_neighbor, vertex])
        right_line = geometry.LineString([right_neighbor, vertex])

        new_left_point = self.get_random_point_on_line(left_line)
        new_right_point = self.find_other_point(new_left_point, left_line, right_line, vertex)

        new_exterior = list(stone.polygon.exterior.coords)
        new_exterior[index] = new_right_point
        new_exterior.insert(index, new_left_point)

        polygon = geometry.Polygon(new_exterior)
        return Stone(polygon)

    def find_other_point(new_left_point, left_line, right_line, center_vertex):
        raise "Unimplemented"

    def get_random_point_on_line(self, line):
        return line.interpolate(self.get_random_number(0.0, 1.0), True).coords[0]

    def generate_random_number(self):
        return random.normalvariate(self.fraction_mean, self.fraction_std)

    def random_vertex(self, polygon, index = None):
        angle = random.uniform(0.0, 2*math.pi)
        center = polygon.centroid.coords[0]
        line = self.create_line_from_center(angle, center)

        while not line.crosses(vertex_line):
            line = affinity.scale(line, xfact=2.0, yfact=2.0, origin=center)

        for i in xrange(len(polygon.exterior.coords)-1):
            coord1 = polygon.exterior.coords[i]
            coord2 = polygon.exterior.coords[i+1]

            polygon_line = geometry.LineString([coord1, coord2])
            if polygon_line.intersects(line):
                intersection = geometry.Point(polygon_line.intersection(line).coords[0])
                if intersection.distance(geometry.Point(coord1)) < intersection.distance(geometry.Point(coord2)):
                    return (polygon.exterior.coords[i], i)
                else:
                    return (polygon.exterior.coords[i+1], i+1)

# This vertex chipper selects a random vertex to chip, then gets the lines that
# form that vertex. It picks one point from each intersecting line and connects
# those two points to chip the stone.
#
# The points are chosen with a normal distribution centered about the
# +fraction_mean+ of the line segment. The standard deviation is given by the 
# +fraction_std+.
#
# A +fraction_mean+ of 0.5 means that the normal distribution is centered around
# the midpoint of the line segment. A +fraction_mean+ of 0.1 means that the normal
# distribution is centered around 0.1 times the distance of segment.
class RandomPointVertexChipper(VertexChipper):
    def __init__(self, fraction_mean = 0.5, fraction_std = 0.15):
        self.fraction_mean = fraction_mean
        self.fraction_std = fraction_std

    def chip_vertex(self, stone, index = None):
        vertex, index = self.random_vertex(stone.polygon, index)
        left_neighbor = stone.polygon.exterior.coords[index-1]
        right_neighbor = stone.polygon.exterior.coords[index+1]
        left_line = geometry.LineString([left_neighbor, vertex])
        right_line = geometry.LineString([right_neighbor, vertex])

        new_left_point = self.get_random_point_on_line(left_line)
        new_right_point = self.get_random_point_on_line(right_line)

        new_exterior = list(stone.polygon.exterior.coords)
        new_exterior[index] = new_right_point
        new_exterior.insert(index, new_left_point)

        polygon = geometry.Polygon(new_exterior)
        return Stone(polygon)

    def get_random_point_on_line(self, line):
        return line.interpolate(self.get_random_number(0.0, 1.0), True).coords[0]

    def generate_random_number(self):
        return random.normalvariate(self.fraction_mean, self.fraction_std)

# This vertex chipper creates chips by selecting random angles. The starting angle
# is chosen uniformly from the range [0, 2*pi]. Then a second angle is chosen from
# the range [0, pi] from a normal distribution centered around mean +angle_mean+
# with standard deviation of +angle_std+.
#
# We now have two angles, starting_angle and (starting_angle + second_angle). These
# two angles are defined with respect to the centroid of the stone. A straight
# line is drawn between the intersection of the first and second angles with the
# stone. Vertices that lie outside of this new line are removed from the stone.
class AngleVertexChipper(VertexChipper):
    def __init__(self, angle_mean, angle_std):
        self.angle_mean = angle_mean
        self.angle_std = angle_std

    def chip_vertex(self, stone, index = None):
        centroid = stone.polygon.centroid.coords[0]
        angle = self.get_random_starting_angle()
        first_line = self.create_line_from_center(angle, centroid)

        angle_change, is_clockwise = self.get_random_angle()
        second_line = self.create_rotated_line(first_line, angle_change)

        polygon = self.remove_vertices_between(stone.polygon, first_line, second_line, is_clockwise)
        return Stone(polygon)

    # Removes the vertices between +point+ and +vertex+.
    def remove_vertices_between(self, polygon, first_line, second_line, is_clockwise=True):
        coordinates = list(polygon.exterior.coords)
        if is_clockwise:
            coordinates.reverse()

        first_index = None
        second_index = None
        for i in xrange(len(coordinates)-1):
            current_vertex = coordinates[i]
            next_vertex = coordinates[i+1]
            line = geometry.LineString([current_vertex, next_vertex])

            if line.intersects(first_line):
                first_index = i
                first_point = line.intersection(first_line)
            if line.intersects(second_line):
                second_index = i
                second_point = line.intersection(second_line)

        inserted_points = [first_point.coords[0], second_point.coords[0]]
        new_coordinates = self.create_new_coordinates(coordinates, first_index, second_index, inserted_points)
        if is_clockwise:
            new_coordinates = new_coordinates.reverse()

        return geometry.Polygon(new_coordinates)

    def create_new_coordinates(self, coordinates, first_index, second_index, inserted_points):
        if second_index >= first_index:
            return coordinates[:(first_index+1)] + inserted_points + coordinates[(second_index+1):]
        else:
            return inserted_points + coordinates[(second_index+1):(first_index+1)] + [inserted_points[0]]

    def almost_equals(self, first_tuple, second_tuple):
        first_point = geometry.Point(first_tuple)
        second_point = geometry.Point(second_tuple)
        return first_point.almost_equals(second_point)

    # Returns and angle and whether or not you are moving in the forward direction
    # in the clockwise direction
    def get_random_angle(self):
        angle = self.get_random_number(0.0, math.pi)
        if angle < math.pi:
            return (angle, False)
        else:
            return (angle, True)

    def get_random_starting_angle(self):
        return random.uniform(0.0, 2*math.pi)

    def generate_random_number(self):
        return random.normalvariate(self.angle_mean, self.angle_std)

class Stone:
    def __init__(self, polygon):
        self.polygon = polygon
