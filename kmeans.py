import math

class Point:
    def __init__(self, data) -> None:
        self.data = data
        self.cluster_index = -1

    def __str__(self) -> str:
        return self.data.__str__()

    def distance(self, other):
        sum = 0
        for j in range(len(self.data)):
            sum += (self.data[j] - other.data[j]) ** 2
        return math.sqrt(sum)

    def match_cluster(self, list_of_clusters):
        min_distance = self.distance(list_of_clusters[0].centroid)
        nearest_cluster = list_of_clusters[0]
        for cur_cluster in list_of_clusters:
            distance = self.distance(cur_cluster.centroid)
            if distance < min_distance:
                min_distance = distance
                nearest_cluster = cur_cluster
        return nearest_cluster


class Cluster:
    def __init__(self, centroid, index) -> None:
        self.centroid = centroid
        self.points = []
        self.index = index

    def __str__(self) -> str:
        return self.centroid.__str__()

    def add_point(self, point):
        self.points.append(point)

    def clear(self):
        self.points.clear()

    def update_centroid(self, eps):
        mean = [0] * len(self.centroid.data)
        for i in range(len(self.centroid.data)):
            for p in self.points:
                mean[i] += p.data[i]
            mean[i] /= len(self.points)
        new_centroid = Point(mean)
        converged = new_centroid.distance(self.centroid) < eps
        self.centroid = new_centroid
        return converged

def fit(datapoints, k, eps, iter):
    points = [Point(data) for data in datapoints]
    clusters = []
    for i in range(k):
        clusters.append(Cluster(points[i], i))
    i = 0
    converged = []
    for i in range(iter):
        for point in points:  # adding each point to the nearest cluster's points list
            nearest_cluster = point.match_cluster(clusters)
            nearest_cluster.add_point(point)
            point.cluster_index = nearest_cluster.index
        for cluster in clusters:
            converged.append(cluster.update_centroid(eps))
            cluster.points.clear()
        if False not in converged:
            break
        else:
            converged.clear()
    """
    for cluster in clusters:
        cluster.centroid.data = ["%.4f" % i for i in cluster.centroid.data]
        print(*cluster.centroid.data, sep=",")
    """
    return [p.cluster_index for p in points]