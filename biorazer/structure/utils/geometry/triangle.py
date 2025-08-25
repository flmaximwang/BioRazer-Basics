import numpy as np


def edge_lengths(p1, p2, p3):
    a = np.linalg.norm(p2 - p3)
    b = np.linalg.norm(p1 - p3)
    c = np.linalg.norm(p1 - p2)
    return a, b, c


def circumcenter(p1, p2, p3):
    a, b, c = edge_lengths(p1, p2, p3)
    a2 = a * a
    b2 = b * b
    c2 = c * c
    a4 = a2 * a2
    b4 = b2 * b2
    c4 = c2 * c2
    p = a2 * (b2 + c2 - a2) * p1 + b2 * (c2 + a2 - b2) * p2 + c2 * (a2 + b2 - c2) * p3
    p /= 2 * (a2 * b2 + b2 * c2 + c2 * a2) - a4 - b4 - c4
    return p


def circumradius(p1, p2, p3):
    a, b, c = edge_lengths(p1, p2, p3)
    area = 0.5 * np.linalg.norm(np.cross(p2 - p1, p3 - p1))
    return (a * b * c) / (4 * area)


def incenter(p1, p2, p3):
    a, b, c = edge_lengths(p1, p2, p3)
    p = (a * p1 + b * p2 + c * p3) / (a + b + c)
    return p
