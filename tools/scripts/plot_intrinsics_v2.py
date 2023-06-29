import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse

font = {'family' : 'normal',
        # 'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)


parser = argparse.ArgumentParser(description="Generate images from perf tests output")
parser.add_argument("--testteefile", type=str, default="quality/runoutputs/test_basic_perf.txt")
parser.add_argument("--testsgxfile", type=str, default="quality/runoutputs/test_sgx_perf.txt")
parser.add_argument("--outfile", type=str, default="quality/graphs/intrinsics.png")

args = parser.parse_args()

# UNDONE(): confirm enc and dec (i think they should be smaller).
#
# names:
# mov, disk, enc, dec, ocall, ewb
#

iops = 220000
disk = 1.0 / iops *  10**9
# This data is for Intel Xeon E-2200
#
data_v2 = (232, disk, 1599, 1736, 7565, 12265)

def make_pts(mov, disk, enc, dec, ocall, ewb):
    # stack: enc/mov, dec/mov, ocall/mov
    return [1, enc/mov, dec/mov, disk/mov, ocall/mov, ewb/mov, ocall/mov, ewb/mov]

pts_v1 = make_pts(*data_v2)
xaxis = ["Mov", "Enc", "Dec", "Disk", "PageSwap-OCall", "PageSwap-EWB", "PageSwap-OCall-Disk", "PageSwap-EWB-Disk"]
fig, ax = plt.subplots()


p = pts_v1[:len(xaxis)]
p[1] = p[2] = p[3] = 0
ax.barh(xaxis, p, height=0.8)
pp = p[:]
p = np.zeros(len(xaxis))
p[1] += pts_v1[1]
p[4] += pts_v1[1]
p[6] += pts_v1[1]
ax.barh(xaxis, p, height=0.8, left=pp)
pp += p

p = np.zeros(len(xaxis))
p[2] += pts_v1[2]
p[4] += pts_v1[2]
p[6] += pts_v1[2]
ax.barh(xaxis, p, height=0.8, left=pp)
pp += p

p = np.zeros(len(xaxis))
p[3] += pts_v1[3]
p[6] += pts_v1[3]
p[7] += pts_v1[3]
ax.barh(xaxis, p, height=0.8, left=pp)
pp += p

p = [0.0 for _ in range(len(xaxis))]
p[0]=1
ax.barh(xaxis, p, height=0.8)


ax.set_xlabel("Time relative to Mov")
# ax.set_ylabel("Operation Types")
fig.gca().set_aspect(6.0)
fig.subplots_adjust(left=0.3)
plt.show()
