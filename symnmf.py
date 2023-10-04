import sys
import math
import numpy as np
import symNMF

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

try:
    np.random.seed(0)
    k = int(sys.argv[1])
    goal = sys.argv[2]
    filename = sys.argv[3]
    dp = read_from_file(filename)
    result = dp
    if goal=='symnmf':
        w = symNMF.norm(dp)
        h = init_h(k, len(dp), w)
        result = symNMF.symnmf(h, w)
    elif goal=='sym':
        result = symNMF.sym(dp)
    elif goal=='ddg':
        result = symNMF.ddg(dp)
    elif goal=='norm':
        result = symNMF.norm(dp)
    elif goal=='init':
        w = symNMF.norm(dp)
        result = init_h(k, len(dp), w)
    # print result
    for row in result:
        str_list = ["%.4f" % cell for cell in row]
        print(*str_list, sep=',')
except:
    handle_error()
