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
,10,100,1000
]
# ,8000]
# pathoram, ringoram, pathoram_enc = odsl_parse_file(args.testteefile)
# pathoram = odsl_parse_file(args.testteefile)
# signal = [signal_parse_file(args.signalteefile + str(i) + ".txt") for i in signal_query_sizes]


y_pathoram_l1 = [37884,49444,55439,70629,81214,123460,111084,181579,166483,177772,191968,346514,388553,558952,571109,1094037,1754937,3059399,4004581]
y_pathoram_l2 = [37293, 39283,50527,56818,71903,83086,103325,139273,144149,150481,195342,201938,218551,318018,330412,859201,1110275,1670617,2994944]
y_pathoram_l3 = [37326,51085,55751,72568,80910,104677,111526,138365,148118,176290,191352,217891,231935,324047,321297,605285,1088383,1260001,1633021]
y_pathoram_l4 = [32420,43282,48803,58952,71450,90898,98222,117762,130537,165392,167661,187979,205651,230658,285333,655870,1039932,1021986,1107804]
y_pathoram_l5 = [31424,37840,50310,57426,70176,81531,103467,178283,171035,163639,180765,194289,205295,363492,439058,583560,1182627,1797427,2319569]
pts_pathoram_l1 = [(256<<i, y_pathoram_l1[i]/1000.0) for i in range(len(y_pathoram_l1))]
pts_pathoram_l2 = [(256<<i, y_pathoram_l2[i]/1000.0) for i in range(len(y_pathoram_l2))]
pts_pathoram_l3 = [(256<<i, y_pathoram_l3[i]/1000.0) for i in range(len(y_pathoram_l3))]
pts_pathoram_l4 = [(256<<i, y_pathoram_l4[i]/1000.0) for i in range(len(y_pathoram_l4))]
pts_pathoram_l5 = [(256<<i, y_pathoram_l5[i]/1000.0) for i in range(len(y_pathoram_l5))]
colors1 = iter(["aqua","dodgerblue","blue","darkblue"][::-1]) #,"midnightblue"
colors2 = iter(["coral", "red", "firebrick", "maroon", "purple"][::-1]) #, "brown"

# for i in range(len(signal_query_sizes))[::-1]:
  # qs = signal_query_sizes[i]
  # print(signal[i])
  # plt.plot(*zip(*signal[i]), "-^", color=next(colors1), label=f"Signal SGXv2 $\\beta$={qs}", linewidth=1)  

# If we want to show the plots at the multi query levels
def m(l, mult):
    return [(x[0], x[1]*mult) for x in l]


plt.plot(*zip(*pts_pathoram_l1), "-o", color=next(colors2), label=f'page=296B', linewidth=1)
plt.plot(*zip(*pts_pathoram_l2), "-o", color=next(colors2), label=f'page=824B', linewidth=1)
plt.plot(*zip(*pts_pathoram_l3), "-o", color=next(colors2), label=f'page=1880B', linewidth=1)
plt.plot(*zip(*pts_pathoram_l4), "-o", color=next(colors2), label=f'page=3992B', linewidth=1)
plt.plot(*zip(*pts_pathoram_l5), "-o", color=next(colors2), label=f'page=8216B', linewidth=1)

# for i in range(len(signal_query_sizes))[::-1]:
#     pts_pathoram_this = m(pts_pathoram_v1,qs)
#     plt.plot(*zip(*pts_pathoram_this), "-o", color='green', label=f'ENIGMAP SGXv1 $\\beta$={qs}', linewidth=1)



  
plt.xscale("log", base=2)
# plt.yscale("log", base=2)
# plt.yticks([(1<<i) for i in range(0,26,2)])
# plt.xticks([10**i for i in range(1,9,1)])
# plt.xticks([5*10**6, 10**7, 5*(10**7,) 10**8])
plt.tick_params(axis='x', which='minor')
locmin = matplotlib.ticker.LogLocator(base=2.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
plt.gca().get_xaxis().set_minor_locator(locmin)
plt.gca().get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
plt.xlim(2**19.99, 1*10**8)
plt.ylim(2**7, 2**12)
plt.xlabel("Database size")
plt.ylabel("Query Time ($\mu s$)")
plt.axvline(x=2**21, color='gray', ls=':')
plt.text(1.03*2**21,1300,'RAM Swap',rotation=90)
plt.axvline(x=2**23, color='black', ls=':')
plt.text(1.03*2**23,3000,'Disk Swap',rotation=90)

# plt.legend(loc="best",bbox_to_anchor =(1.0,0.85))
plt.legend(loc=2)
# plt.show()
plt.savefig(args.outfile, bbox_inches='tight')