import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse

font = {'family' : 'normal',
        # 'weight' : 'bold',
        'size'   : 14}

matplotlib.rc('font', **font)

parser = argparse.ArgumentParser(description="Generate images from perf tests output")
parser.add_argument("--testteefile", type=str, default="quality/runoutputs/test_basic_perf.txt")
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
    ret += [(n, t)]
  return ret

signal_query_sizes = [1]
# ,8000]
# pathoram, ringoram, pathoram_enc = odsl_parse_file(args.testteefile)
# pathoram = odsl_parse_file(args.testteefile)
signal = [signal_parse_file(args.signalteefile + str(i) + ".txt") for i in signal_query_sizes]
signal_y = [3.5102,5.7799,10.447,20.510,41.824,81.501,174.19,295.81,728.02,1500,3130,5740,12920,26750,51310,86850,174090,347020,687180,1368400,2759700,5516048,10867532,218337330,448260220]
signal_data = [(256<<i, signal_y[i]) for i in range(len(signal_y))]

# pts_pathoram = [(x["oram.N"], x["us/query"]) for x in pathoram]
y_pathoram = [33762,45776,56229,65422,69761,93655,105068,147347,151069,159308,166777,232594,242604,255246,263556,338857,346972,362525,652069,1068319,1153113,1337536,1469272,1805145,1898574]

pts_pathoram = [(256<<i, y_pathoram[i]/1000.0) for i in range(len(y_pathoram))]

y_signalv2 = [133749,98257,601450,543808,859014,608530,962152,815432,1418503,1312416,1838202,2158955,2772744]

# UNDONE: These take too long to initialize, the following point is got by not initiliazing and hardcoding the cycles constant based on the previous rounds.
#
y_signalv2 += [7211675, 9397565, 12320295, 16307730, 22565235, 27438447, 38408868]

pts_signalv2 = [(2<<(i+14), y_signalv2[i]/1000.0) for i in range(len(y_signalv2))]


colors1 = iter(["aqua","dodgerblue","blue","darkblue"][::-1]) #,"midnightblue"
colors2 = iter(["coral", "red", "firebrick", "maroon"][::-1]) #, "brown"

# Predict extra points for signalv1:
#
for i in range(len(signal_query_sizes))[::-1]:
  signal[i] += [(signal[i][-1][0]*2, signal[i][-1][1]*2)]
  signal[i] += [(signal[i][-1][0]*2, signal[i][-1][1]*2)]
for i in range(len(signal)):
  base = 31 # 28 or 32
  for j in range(base, len(signal[i])):
    signal[i][j] = (signal[i][j][0], signal_data[j-8][1])

plt.plot(*zip(*signal[i]), "-^", color='aqua', label=f"Signal", linewidth=1)
plt.plot(*zip(*pts_signalv2), "-^", color='green', label=f"Signal HT (2022)", linewidth=1)
plt.plot(*zip(*pts_pathoram), "-o", color='coral', label=f'ENIGMAP', linewidth=1)
    # if not args.showguidelines:
    #   break


  
plt.xscale("log", base=10)
plt.yscale("log", base=2)
plt.yticks([(1<<i) for i in range(0,29,2)])
plt.xticks([10**i for i in range(1,10,1)])
plt.tick_params(axis='x', which='minor')
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
plt.gca().get_xaxis().set_minor_locator(locmin)
plt.gca().get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())

plt.xlim(1.2*2**14, 4.5*10**9)
plt.ylim(4, 2**29)


# SGXv2:
plt.axvline(x=2**26, color='gray', ls=':')
plt.text(1.03*2**26,4*1.3,'RAM Swap',rotation=90)
plt.axvline(x=2**27, color='black', ls=':')
plt.text(1.03*2**27,4*1.3,'Disk Swap',rotation=90)

plt.axvline(x=2**30, color='gray', ls='-.')
plt.text(1.03*2**30,4*1.3,'Signal RAM Swap',rotation=90)
plt.axvline(x=2**31, color='black', ls='-.')
plt.text(1.03*2**31,4*1.3,'Signal Disk Swap',rotation=90)

plt.axvline(x=2**28, color='gray', ls='-.')
plt.text(1.03*2**28,4*1.3,'Signal HT RAM+Disk Swap',rotation=90)


plt.xlabel("Database size")
plt.ylabel("Query Time ($\mu s$)")
plt.legend(loc="best",bbox_to_anchor =(1.0,0.85))
# plt.show()
plt.savefig(args.outfile, bbox_inches='tight')