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

# dataset_insertion = odsl_extract_points_query(args.testinsertionfile)
# dataset_search = odsl_extract_points_query(args.testsearchfile)

# pts_search = [(x["max nodes"], x["us/query"]) for x in dataset_search]
# pts_insertion = [(x["max nodes"], x["us/query"]) for x in dataset_insertion]
# pts_deletion = [(x["max nodes"], 5*x["us/query"]) for x in dataset_insertion]
# search_t = [34712,45748,51406,64286,74915,93221,101029,149997,148202,163164,176706,237723,233980,250022,265588,349174,346239,366388,384583,633583,700413,857737, 961238, 991203, 1045923]
search_t = [33762,45776,56229,65422,69761,93655,105068,147347,151069,159308,166777,232594,242604,255246,263556,338857,346972,362525,652069,1068319,1153113,1337536,1469272,1805145,1898574]
insertion_t = [69423,91495,102811,128572,149829,186442,202057,299994,296404,326328,353411,475446,467959,500044,531175,698348,692478,732776,1307580,2154182,2381402,2916304,3268207,3370088,2091845]

pts_search = [(256<<i, search_t[i]/1000) for i in range(len(search_t))]
pts_insertion = [(256<<i, insertion_t[i]/1000) for i in range(len(insertion_t))]
pts_deletion = [(256<<i, 5*insertion_t[i]/1000) for i in range(len(insertion_t))]

plt.plot(*zip(*pts_deletion), "-o", color="red", label=f'Deletion', linewidth=1)
plt.plot(*zip(*pts_insertion), "-o", color="purple", label=f'Insertion', linewidth=1)
plt.plot(*zip(*pts_search), "-o", color="blue", label=f'Search', linewidth=1)

plt.xscale("log", base=10)
plt.yscale("log", base=10)
plt.xlim(2**17, 2**31)
plt.ylim(0.9*10, 2*10**4)
plt.xlabel("Database size")
plt.ylabel("Query time($\mu s$)")
plt.legend(loc="best")
plt.axvline(x=2**26, color='gray', ls=':')
plt.text(1.03*2**26,2**4,'RAM Swap',rotation=90)
plt.axvline(x=2**27, color='black', ls=':')
plt.text(1.03*2**27,2**4,'Disk Swap',rotation=90)
# plt.show()
plt.savefig(args.outfile, bbox_inches='tight')
