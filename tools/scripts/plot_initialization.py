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
parser.add_argument("--testinitfile", type=str, default="initialization_sample.txt")
parser.add_argument("--testinsfile", type=str, default="insert_sample.txt")
parser.add_argument("--testbitinitfile", type=str, default="initialization_bitonic.txt")
parser.add_argument("--outfile", type=str, default="quality/graphs/initialization.png")

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

def add_external_memory_cost_to_bintonic(pts):
  return [
    (
      x[0], 
      x[1]*(1 + 12.5*(max(0,x[0]-(2**20))/x[0]))
    )
    for x in pts
  ]

dataset = odsl_extract_points_query(args.testinitfile)
dataset_bitonic = odsl_extract_points_query(args.testbitinitfile)
dataset_insertion = odsl_extract_points_query(args.testinsfile)

pts = [(x["max nodes"], x["s initialization time"]) for x in dataset]
pts_insertion = [(x["max nodes"], x["max nodes"]*x["us/query"]/1000000.0) for x in dataset_insertion]
pts_bitonic = [(x["max nodes"], x["s initialization time"]) for x in dataset_bitonic]
pts_bitonic = add_external_memory_cost_to_bintonic(pts_bitonic)

plt.plot(*zip(*pts_bitonic), "-o", color="red", label=f'Oblix Initialization', linewidth=1)
plt.plot(*zip(*pts_insertion), "-o", color="purple", label=f'Naive Initialization', linewidth=1)
plt.plot(*zip(*pts), "-o", color="blue", label=f'Fast Initialization', linewidth=1)

plt.xscale("log", base=10)  
plt.yscale("log", base=10)
plt.xlim(10**2, 2*2**28)
plt.ylim(0.01, 7*10**5)
plt.xlabel("Database size")
plt.ylabel("Initialization Time (s)")
plt.legend(loc="best")#,bbox_to_anchor =(1.0,0.65))
plt.show()
plt.savefig(args.outfile, bbox_inches='tight')
