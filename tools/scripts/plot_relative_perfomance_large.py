import matplotlib
import matplotlib.pyplot as plt

font = {'family' : 'normal',
        # 'weight' : 'bold',
        'size'   : 14}
matplotlib.rc('font', **font)

width = 0.35
labels = ['Base', '+Locality-Friendly\n Layout', '+Page\nCache', '+Bucket\nCache']
overhead_types = ['Computation', 'OCall+Enc+Dec', 'Disk IO']

# oram depth:
#
depth = 20
callsperaccess = 1
# 0.984 comes from the profiling_test.cpp, it is the cost to ocall once.
# We do callsperaccess ocalls per oram access
# (depends on the oram depth it is ceil((Z*blocksize*logN)/pagesize).
#
ocalltime = depth*callsperaccess*0.984



# Actual performance number of one query with all optimizations, with profiling disabled (us)
#

def get_data(t, n, l, cl0, cl1, fc):
  nn = 1.44 * n
  reads = nn
  ocalls = (reads-cl0) * (n - cl1) / l / 2
  to = 7.5 * ocalls
  tio = ocalls * (1000000/220000) * ((n-fc) / n) * 2
  tins = t - to - tio
  print(n, nn, reads, ocalls)
  return tins, to, tio


# For large databases the numbers come from the time spent inside the enclave vs the time spent in
# ocalls. This is obtained by counting the number of ocalls, and comparing the time 
# of a version that uses in memory server with a version that uses disk.
#
n = 28
l = 4

# No io points:
#
points = [
  get_data(2618.5, 26, 1.5, 4, 0, 26),
  get_data(1212.061, 26, 4, 4, 0, 26),
  get_data(723.061, 26, 4, 20, 0, 26),
  get_data(542.061, 26, 4, 20, 10, 26),
]

# Always io points:
#
#
points_ = [
  get_data(3118.5, 28, 1.8, 4, 10, 3),
  get_data(1512.061, 28, 4, 4, 10, 3),
  get_data(923.061, 28, 4, 18, 10, 10),
  get_data(605.061, 28, 4, 18, 18, 10),
]
  

print(points)
fig, ax = plt.subplots()
pts_prev = [0 for i in range(len(points))]
for i, tp in enumerate(overhead_types):
  pts = [points[j][i]/(points[-1][0]+points[-1][1]+points[-1][2]) for j in range(len(points))]
  ax.bar(labels, pts, width, label=tp, bottom=pts_prev)
  pts_prev = [pts_prev[j] + pts[j] for j in range(len(points))]
ax.set_ylabel('Relative overhead')
# ax.set_title('Optimization impact analysis for $n = 2^{20}$')
ax.legend()

ax.set_ybound(0,6.0)

plt.show()