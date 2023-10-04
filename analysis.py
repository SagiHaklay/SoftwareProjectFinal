import sys
import sklearn.metrics as metrics
import kmeans
import symNMF
import numpy as np
import math

def read_from_file(filename):
    datapoints = []
    f = open(filename, "r")
    for line in f:
        curr_line = line.split(',')
        curr_line = [float(x) for x in curr_line]
        datapoints.append(curr_line)
    return datapoints

def handle_error():
    print('An error has occurred')
    sys.exit()

def init_h(k, n, w):
    m = np.average(w)
    h = [[np.random.uniform(low=0, high=2*math.sqrt(m/k)) for j in range(k)] for i in range(n)]
    return h

def to_tag_array(h):
    tags = []
    for row in h:
        maxi = 0
        for i in range(1, len(row)):
            if row[i] > row[maxi]:
                maxi = i
        tags.append(maxi)
    return tags

try:
    np.random.seed(0)
    k = int(sys.argv[1])
    filename = sys.argv[2]
    dp = read_from_file(filename)
    w = symNMF.norm(dp)
    h = init_h(k, len(dp), w)
    nmf_result = symNMF.symnmf(h, w)
    kmeans_result = kmeans.fit(dp, k, 0.001, 600)
    nmf_score = metrics.silhouette_score(dp, to_tag_array(nmf_result))
    kmeans_score = metrics.silhouette_score(dp, kmeans_result)
    print('nmf:', "%.4f" % nmf_score)
    print('kmeans:', "%.4f" % kmeans_score)
except:
    handle_error()