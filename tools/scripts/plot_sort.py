import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse
import math

font = {'family' : 'normal',
        # 'weight' : 'bold',
        'size'   : 14}

matplotlib.rc('font', **font)

parser = argparse.ArgumentParser(description="Generate images from perf tests output")
parser.add_argument("--infile", type=str, default="quality/runoutputs/test_basic_perf.txt")
parser.add_argument("--outfilepages", type=str, default="quality/graphs/sort_pages.png")
parser.add_argument("--outfileinstr", type=str, default="quality/graphs/sort_instr.png")
parser.add_argument("--outfiletime", type=str, default="quality/graphs/sort_time.png")

args = parser.parse_args()

def parse_report(report: str):
  res = {}
  for line in report.split("\n"):
    if line.count(":") == 1:
      k, v = line.split(":")
      v = float(v.strip())
      res[k] = v
  return res


def extract_points_from_file(input_file: str):
  ret1, ret2 = [], []
  with open(input_file, "r") as f:
    s = f.read()
    s = s.split("] TestSort.TestBucketObliviousSortPerf")[1]
  reports = s.split("test bucket oblivious sort perf")[1:]

  for report in reports:
    md, reps = report.split("Running bucket oblivious sort")
    r1, r2 = reps.split("Running bitonic sort")
    sz = int(md.split("\n")[0].strip())
    r1 = r1.strip()
    r2 = r2.strip()
    pts1 = parse_report(r1)
    pts2 = parse_report(r2)
    
    ret1 += [pts1]
    ret2 += [pts2]
    
  return ret1, ret2

bos, bits = extract_points_from_file(args.infile)

pts_time = [(x["N"], x["time (s)"]) for x in bos]
pts_access_count = [(x["N"], x["Access count"]) for x in bos]
pts_swap_count = [(x["N"], x["Swap count"]) for x in bos]
pts_page_swap = [(x["N"], x["Write pages"]) for x in bos]

plt.plot(*zip(*pts_page_swap), "-o", color="red", label=f'Bucket Oblivious Sort', linewidth=1)
plt.xscale("log", base=10)
plt.yscale("log", base=10)
# plt.xlim(10**2, 3*2**26)
# plt.ylim(0.9*10, 2*10**4)
plt.xlabel("Vector size")
plt.ylabel("Sorting page swaps")
plt.legend(loc="best")
# plt.show()
plt.savefig(args.outfilepages, bbox_inches='tight')
plt.close()

plt.plot(*zip(*pts_access_count), "-o", color="red", label=f'Bucket Oblivious Sort Entry Reads', linewidth=1)
plt.plot(*zip(*pts_swap_count), "-o", color="red", label=f'Bucket Oblivious Sort Entry Swaps', linewidth=1)
plt.xscale("log", base=10)
plt.yscale("log", base=10)
# plt.xlim(10**2, 3*2**26)
# plt.ylim(0.9*10, 2*10**4)
plt.xlabel("Vector size")
plt.ylabel("Count")
plt.legend(loc="best")
# plt.show()
plt.savefig(args.outfileinstr, bbox_inches='tight')
plt.close()

plt.plot(*zip(*pts_time), "-o", color="red", label=f'Bucket Oblivious Sort Time', linewidth=1)
plt.xscale("log", base=10)
plt.yscale("log", base=10)
# plt.xlim(10**2, 3*2**26)
# plt.ylim(0.9*10, 2*10**4)
plt.xlabel("Vector size")
plt.ylabel("Time (s)")
plt.legend(loc="best")
# plt.show()
plt.savefig(args.outfiletime, bbox_inches='tight')
plt.close()