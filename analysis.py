import sys
import sklearn.metrics as metrics
import kmeans


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

k = int(sys.argv[1])
filename = sys.argv[2]
dp = read_from_file(filename)
kmeans_result = kmeans.fit(dp, k, 0.001, 600)
print(kmeans_result)
kmeans_score = metrics.silhouette_score(dp, kmeans_result)
print(kmeans_score)
# Calculate H with SymNMF, create label array, and get silhouette score