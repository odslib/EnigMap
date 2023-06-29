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

# dataset = odsl_extract_points_query(args.testinitfile)
# dataset_bitonic = odsl_extract_points_query(args.testbitinitfile)
# dataset_insertion = odsl_extract_points_query(args.testinsfile)

# pts = [(x["max nodes"], x["s initialization time"]) for x in dataset]
# pts_insertion = [(x["max nodes"], x["max nodes"]*x["us/query"]/1000000.0) for x in dataset_insertion]
# pts_bitonic = [(x["max nodes"], x["s initialization time"]) for x in dataset_bitonic]
# pts_bitonic = add_external_memory_cost_to_bintonic(pts_bitonic)


insertion_t = [69423,91495,102811,128572,149829,186442,202057,299994,296404,326328,353411,475446,467959,500044,531175,698348,692478,732776,1307580,2154182,2381402,2916304,3268207,3370088,2091845]

y = [213619897,302907758,106690329,630547621,187154011,237458975,1134811856,2469465211,4343608189,7930971719,25388398818,68790097216,107168702034,290801719093,655047632152,979591415786,2851959260903,8548245005010,34260169633454,77291307332391,182671752943826,436316718275895]
y_bitonic=[6021120,17031168,46448640,123002880,317915136,805109760,2003828736,4913233920,11890851840,28449964032,67381493760,158148329472,368176005120,850856509440,1953270595584,4456867430400,10113037369344,188304355328000,612980156416000]


pts = [((1<<(8+i)), y[i] / 10**9) for i in range(len(y))]
pts_insertion = [((1<<(8+i)), (1<<(8+i))*insertion_t[i] / 10**9) for i in range(len(insertion_t))]
pts_bitonic = [((1<<(8+i)), y_bitonic[i] / 10**9) for i in range(len(y_bitonic))]

print(pts_insertion)

plt.plot(*zip(*pts_bitonic), "-o", color="red", label=f'Oblix Initialization', linewidth=1)
plt.plot(*zip(*pts_insertion), "-o", color="purple", label=f'Naive Initialization', linewidth=1)
plt.plot(*zip(*pts), "-o", color="blue", label=f'Fast Initialization', linewidth=1)

plt.axvline(x=2**24, color='gray', ls=':')
plt.text(1.03*2**24,1.3,'RAM Swap',rotation=90)
plt.axvline(x=2**25, color='black', ls=':')
plt.text(1.03*2**25,1.3,'Disk Swap',rotation=90)

plt.xscale("log", base=10)  
plt.yscale("log", base=10)
plt.xlim(10**4, 2*2**29)
plt.ylim(0.01, 7*10**6)
plt.xlabel("Database size")
plt.ylabel("Initialization Time (s)")
plt.legend(loc="best")#,bbox_to_anchor =(1.0,0.65))
# plt.show()
plt.savefig(args.outfile, bbox_inches='tight')
