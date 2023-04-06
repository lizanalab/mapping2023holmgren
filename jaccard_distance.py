import numpy as np


def jaccard_distance(c1, c2):
    union = np.union1d(c1, c2).size
    if union < 1:
        return 1
    return 1 - np.intersect1d(c1, c2).size / union


def distance_metric(p1, p2):
    """
    p1 = [
        [1, 2, 3], # community 1
        [4, 5, 6], # community 2
        # ...
        [10, 11] # community n

    ]

    p2 = [...] # as p1
    """

    p1 = np.array(p1)
    p2 = np.array(p2)
    p1_size = p1.size
    p2_size = p2.size

    def dist_impl(p1, p2):
        total = 0

        for c1 in p1:
            c1 = np.array(c1)
            min_dist = np.inf

            for c2 in p2:
                dist = jaccard_distance(c1, c2) * c1.size / p1_size
                min_dist = min(min_dist, dist)

            total += min_dist

        return total

    return 0.5 * dist_impl(p1, p2) + 0.5 * dist_impl(p2, p1)


def test_distance():
    p1 = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
    ]

    p2 = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8],
        [9]
    ]

    return distance_metric(p1, p1), distance_metric(p2, p2), distance_metric(p1, p2), distance_metric(
        [[i] for i in range(10)], [list(range(10))])
