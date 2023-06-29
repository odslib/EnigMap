import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse

font = {'family' : 'normal',
        # 'weight' : 'bold',
        'size'   : 14}

matplotlib.rc('font', **font)

parser = argparse.ArgumentParser(description="Generate images from perf tests output")
parser.add_argument("--testteefile", type=str, default="tools/scripts/data2/infiniteepc.txt")
parser.add_argument("--outfile", type=str, default="quality/graphs/perf.png")
parser.add_argument("--signalteefile", type=str, default="tools/scripts/signal_data/data_")
parser.add_argument("--showguidelines", default=False, action='store_true')

args = parser.parse_args()

def odsl_parse_file(input_file):
  with open(input_file, "r") as f:
    data = f.read()  
  sets = data.split("[----------] 4 tests from PerfOTree/")[1:]
  # sets = data.split("[----------] 1 test from PerfOTree/")[1:]
  assert sets[0].startswith("0, where TypeParam = TestParameter<_ORAM::PathORAM")
  assert sets[1].startswith("0 (")
  # assert sets[2].startswith("1, where TypeParam = TestParameter<_ORAM::RingORAM")
  # assert sets[3].startswith("1 (")
  # assert sets[4].startswith("2, where TypeParam = TestParameter<_ORAM::PathORAM")
  # assert sets[5].startswith("2 (")
  pathoram = odsl_extract_points_query(sets[0],0)
  # ringoram = odsl_extract_points_query(sets[2],1)
  # pathoram_enc = odsl_extract_points_query(sets[4],2)
  return pathoram #, ringoram, pathoram_enc

def odsl_extract_points_query(s: str, tid: int):
  s = s.split("[ RUN      ]")[3].strip()
  ret = []
  assert s.startswith(f"PerfOTree/{tid}.PointSearchPartial")
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

def signal_parse_file(input_file):
  with open(input_file, "r") as f:
    data = f.read()
  sets = data.split("\nhash_lookup_one_query_phone_varying_phone_db/")[1:]
  ret = []
  for pt in sets:
    n = int(pt[:30].split(" ")[0])
    t = float(pt.split("\n")[1].split("[")[1].split("]")[0].split(" ")[2])
    tmult = pt.split("\n")[1].split("[")[1].split("]")[0].split(" ")[3]
    if tmult == "s":
      tmult = 1000000
    elif tmult == "ms":
      tmult = 1000
    elif tmult == "us":
      tmult = 1
    else:
      assert False, f"Invalid tmult: {tmult}"
    t = t * tmult
    if n < 200: continue
    ret += [(n, t)]
  return ret

signal_query_sizes = [1
#,10,100,1000
]
# ,8000]
# pathoram, ringoram, pathoram_enc = odsl_parse_file(args.testteefile)
# pathoram = odsl_parse_file(args.testteefile)
signal = [signal_parse_file(args.signalteefile + str(i) + ".txt") for i in signal_query_sizes]


y_pathoram_v1 = [18.0218,22.4032,30.1097,34.734,50.4785,74.1553,93.6424,139.587,157.301,169.651,220.027,234.004,265.761,305.779,341.667,412.003,439.887,437.374,483.075]
y_pathoram_v2 = [37694,50016,56216,70180,81301,101053,132556,144491,150915,197660,203434,218378,261914,276038,289690,352941,361265,387264,407911]
pts_pathoram_v2 = [(256<<i, y_pathoram_v2[i]/1000.0) for i in range(len(y_pathoram_v2))]
pts_pathoram_v1 = [(256<<i, y_pathoram_v1[i]) for i in range(len(y_pathoram_v1))]
colors1 = iter(["aqua","dodgerblue","blue","darkblue"][::-1]) #,"midnightblue"
colors2 = iter(["coral", "red", "firebrick", "maroon"][::-1]) #, "brown"

for i in range(len(signal_query_sizes))[::-1]:
  qs = signal_query_sizes[i]
  print(signal[i])
  plt.plot(*zip(*signal[i]), "-^", color=next(colors1), label=f"Signal SGXv2 $\\beta$={qs}", linewidth=1)  

# If we want to show the plots at the multi query levels
def m(l, mult):
    return [(x[0], x[1]*mult) for x in l]

for i in range(len(signal_query_sizes))[::-1]:
    qs = signal_query_sizes[i]
    pts_pathoram_this = m(pts_pathoram_v2,qs)
    plt.plot(*zip(*pts_pathoram_this), "-o", color=next(colors2), label=f'ENIGMAP SGXv2 $\\beta$={qs}', linewidth=1)
    if not args.showguidelines:
      break

for i in range(len(signal_query_sizes))[::-1]:
    pts_pathoram_this = m(pts_pathoram_v1,qs)
    plt.plot(*zip(*pts_pathoram_this), "-o", color='green', label=f'ENIGMAP SGXv1 $\\beta$={qs}', linewidth=1)



plt.xscale("log", base=10)
plt.yscale("log", base=2)
plt.yticks([(1<<i) for i in range(0,26,2)])
plt.xticks([10**i for i in range(1,9,1)])
plt.tick_params(axis='x', which='minor')
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
plt.gca().get_xaxis().set_minor_locator(locmin)
plt.gca().get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
plt.xlim(100, 4*10**8)
plt.ylim(1, 2**23)
plt.xlabel("Database size")
plt.ylabel("Query Time ($\mu s$)")
plt.legend(loc="best",bbox_to_anchor =(1.0,0.85))
# plt.show()
plt.savefig(args.outfile, bbox_inches='tight')