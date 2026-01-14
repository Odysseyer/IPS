# Copyright (C) 2026 Zhaoqin Xu, Liu Gan, Mingjie Wei
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
from matplotlib.path import Path

class AttrDict(dict): # Convert dictionary key-value pairs to attributes
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self
class DefinedError(Exception):
    pass


 
def Dealine(line: str) -> list: # Process a line of data
    if '#' in line:
        exegesis = line.index('#')
        newline = line[:exegesis].split()
    else:
        newline = line.strip()
        newline = newline.split()
    if len(newline) != 0: # Filter out empty lines
        for x in range(len(newline)):
            newline[x] = str2num(newline[x])
    return newline


def str2num(strings: str): # Convert string to number
    if strings.isalpha():
        return strings
    elif strings.isdigit() or strings.replace('-','').isdigit():
        return int(strings)
    elif '.' in strings and strings.replace('.', '').replace('-', '').isdigit():
        return float(strings)
    else:
        return strings

          
def orientation(p, q, r):
    """Return orientation between three points: 0 - collinear, 1 - clockwise, 2 - counterclockwise"""
    val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1])
    if val == 0:
        return 0
    elif val > 0:
        return 1
    else:
        return 2

def gift_wrapping(points):
    """Find outermost points and form a polygon"""
    if len(points) < 3:
        return list(range(len(points)))

    hull = []
    leftmost_point_index = min(range(len(points)), key=lambda i: points[i][0])
    p_index = leftmost_point_index

    while True:
        hull.append(p_index)
        q_index = 0
        for r_index in range(len(points)):
            if r_index == p_index:
                continue
            o = orientation(points[p_index], points[q_index], points[r_index])
            if o == 2 or (o == 0 and np.linalg.norm(points[r_index] - points[p_index]) > np.linalg.norm(points[q_index] - points[p_index])):
                q_index = r_index
        p_index = q_index
        if p_index == leftmost_point_index:
            break

    return hull

def rotate(points, angle, origin):
    """Rotate point set around origin"""
    angle_rad = np.radians(angle)
    rotation_matrix = np.array([
        [np.cos(angle_rad), -np.sin(angle_rad)],
        [np.sin(angle_rad), np.cos(angle_rad)]
    ])
    return np.dot(points - origin, rotation_matrix) + origin

def translate(points, translation):
    """Translate point set"""
    return points + translation

def is_point_inside_circle(point, circle_center, radius):
    """Check if point is inside circle"""
    return np.linalg.norm(point - circle_center) < radius

def is_polygon_outside_circles(polygon, circles, radius):
    """
    Check if polygon does not intersect with any circles, has no circles inside it,
    and circles with polygon endpoints as centers and radius do not intersect with other circles
    """
    path = Path(polygon)
    for point in polygon:
        for circle_center in circles:
            if is_point_inside_circle(point, circle_center, radius):
                return False
            if path.contains_point(circle_center):
                return False
            if np.linalg.norm(np.array(point) - np.array(circle_center)) < 2 * radius:
                return False
    return True

def is_polygon_within_bounds(polygon, plane_size):
    """Check if polygon is within bounded plane"""
    x_min, y_min = np.min(polygon, axis=0)
    x_max, y_max = np.max(polygon, axis=0)
    return x_min >= 0 and y_min >= 0 and x_max <= plane_size and y_max <= plane_size

def find_valid_placement(hull, points, circles, radius, plane_size, max_attempts=2000):
    """Try to find a valid placement position"""
    origin = np.mean(points[hull], axis=0)
    for _ in range(max_attempts):
        # Random translation
        translation = np.random.uniform(0, plane_size, size=2)
        translated_points = translate(points, translation - origin)
        
        # Random rotation
        angle = np.random.uniform(0, 360)
        rotated_points = rotate(translated_points, angle, translation)

        # Check if intersects with any circles and within bounds
        if circles.size == 0:
            if is_polygon_within_bounds(rotated_points[hull], plane_size):
                return rotated_points
        elif is_polygon_outside_circles(rotated_points[hull], circles, radius) and is_polygon_within_bounds(rotated_points[hull], plane_size):
            return rotated_points

    raise ValueError('No valid placement found')