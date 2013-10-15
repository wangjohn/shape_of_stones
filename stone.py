import polygon
import random

class Stone:
    def __init__(self, polygon, angle_deviation = 15.0):
        self.polygon = polygon
        self.angle_deviation = angle_deviation

    def chip_vertex(self, index = None):
        if not index:
            index = random.randrange(self.polygon.n)
        vertex = self.polygon.get_vertex(index)
        left_neighbor = self.polygon.get_vertex(index-1)
        right_neighbor = self.polygon.get_vertex(index+1)

        max_length = self.max_chip_length(vertex, left_neighbor, right_neighbor)
        angle = self.get_random_angle()

    def generate_new_vertices(self, angle, 
        
    def max_chip_length(self, vertex, left_neighbor, right_neighbor):
        return float(min(vertex.distance(left_neighbor), vertex.distance(right_neighbor))) / 2.0

    # Returns an angle between -90.0 and +90.0 based on a normal distribution
    # centered with mean of 0.0 and standard deviation which is set by an input
    # to the stone.
    def get_random_angle(self):
        angle = self._generate_potential_angle()):
        while angle <= -90.0 or angle >= 90.0:
            angle = self._generate_potential_angle()

        return angle

    def _generate_potential_angle(self):
        return random.normalvariate(0.0, self.angle_deviation)

class Chip:
    def __init__(self, stone):
        self.stone = stone

    def 
