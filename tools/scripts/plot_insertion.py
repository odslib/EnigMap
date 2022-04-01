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
# parser.add_argument("--testteefile", type=str, default="deletions_sample.txt")
parser.add_argument("--testinsertionfile", type=str, default="insertions_sample.txt")
parser.add_argument("--testsearchfile", type=str, default="search_sample.txt")
parser.add_argument("--outfile", type=str, default="quality/graphs/operations.png")

args = parser.parse_args()

def odsl_extract_points_query(input_file: str):
  ret = []
  with open(input_file, "r") as f:
    s = f.read()
  reports = s.split("[Report]")[1:]
  for report in reports:
    data = report.split("[/Report]")[0].strip()
    res = {}
    for line in data.split("\n"):
      if line.count(":") == 1:
        k, v = line.split(":")
        v = float(v)
        res[k] = v
    ret += [res]
  return ret

dataset_insertion = odsl_extract_points_query(args.testinsertionfile)
dataset_search = odsl_extract_points_query(args.testsearchfile)

pts_search = [(x["max nodes"], x["us/query"]) for x in dataset_search]
pts_insertion = [(x["max nodes"], x["us/query"]) for x in dataset_insertion]
pts_deletion = [(x["max nodes"], 5*x["us/query"]) for x in dataset_insertion]

plt.plot(*zip(*pts_deletion), "-o", color="red", label=f'Deletion', linewidth=1)
plt.plot(*zip(*pts_insertion), "-o", color="purple", label=f'Insertion', linewidth=1)
plt.plot(*zip(*pts_search), "-o", color="blue", label=f'Search', linewidth=1)

plt.xscale("log", base=10)
plt.yscale("log", base=10)
plt.xlim(10**2, 3*2**26)
plt.ylim(0.9*10, 2*10**4)
plt.xlabel("Database size")
plt.ylabel("Query time($\mu s$)")
plt.legend(loc="best")
# plt.show()
plt.savefig(args.outfile, bbox_inches='tight')
