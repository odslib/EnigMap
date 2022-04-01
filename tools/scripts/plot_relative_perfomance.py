import matplotlib
import matplotlib.pyplot as plt

font = {'family' : 'normal',
        # 'weight' : 'bold',
        'size'   : 14}
matplotlib.rc('font', **font)

width = 0.35
labels = ['Base', '+Locality-Friendly\n Layout', '+Page\nCache', '+Bucket\nCache']
overhead_types = ['Instruction', 'OCall+Encryption', 'Disk IO']

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
base_divider = 245

# The format is actually: total time, encryption rate, io
# total time is just the time of the get operation
# encryption+ocall is percentage of time spent in encryption
# These numbers come from the profile graphs (improvements_x.cpp)
#
points = [[575,0.56,575-521],[521,0.47,521 - 498],[498,0.45,0],[447,0.40,0]]


# This loop turns the points into the format: [inst time, ocall+enc time, io time]
for i, point in enumerate(points):
  # We need to turn the instruction cost into the actual instruction cost without the profiling overhead
  # to do so, we just subtract the number of events times the profiling cost
  # (alternatively we could divide the instr time by a constant, but it is similar):
  #
  point[0] -= 1000 * 0.2

  # Remove io from total cost:
  #
  point[0] -= point[2]

  # We can't actually measure average ocall time, since we can't measure time
  # inside of enclaves, but it supposed to be pretty small,
  # specifically, ocall/querys is less than one 0.98us, and we should do 20 ocalls per 
  # query, except if packing is enclave cache is enabled.
  #
  point[1] = ocalltime + point[0]*point[1]


  # Remove encryption from total cost:
  #
  point[0] -= point[1]
  

print(points)
fig, ax = plt.subplots()
pts_prev = [0 for i in range(len(points))]
for i, tp in enumerate(overhead_types):
  pts = [points[j][i]/base_divider for j in range(len(points))]
  ax.bar(labels, pts, width, label=tp, bottom=pts_prev)
  pts_prev = [pts_prev[j] + pts[j] for j in range(len(points))]
ax.set_ylabel('Relative overhead')
# ax.set_title('Optimization impact analysis for $n = 2^{20}$')
ax.legend()

ax.set_ybound(0,1.8)

plt.show()