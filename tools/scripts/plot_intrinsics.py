import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse

font = {'family' : 'normal',
        # 'weight' : 'bold',
        'size'   : 16}

matplotlib.rc('font', **font)


parser = argparse.ArgumentParser(description="Generate images from perf tests output")
parser.add_argument("--testteefile", type=str, default="quality/runoutputs/test_basic_perf.txt")
parser.add_argument("--testsgxfile", type=str, default="quality/runoutputs/test_sgx_perf.txt")
parser.add_argument("--outfile", type=str, default="quality/graphs/intrinsics.png")

args = parser.parse_args()

def odsl_parse_file(input_file):
  with open(input_file, "r") as f:
    data = f.read()  
  sets = data.split("[ RUN      ] TT/Perf_Encryption/")[1:]
  sets += data.split("[ RUN      ] TT/BasicPerf/")[1:]
  rets = {}
  for s in sets:
    s = s.split("[       OK ]")[0]
    szlog2 = int(s[:10].split(".")[0])
    testtype = s[:100].split(".")[1].split("\n")[0].split(" ")[0].strip()
    sz = 1 << szlog2
    lat_us = float(s.split("us/op")[0].split("\n")[-1].strip())
    if testtype not in rets:
      rets[testtype] = {}
    rets[testtype][sz] = lat_us
  return rets

def odsl_parse_sgx_files(input_file):
  testtype = 'OCallSwap'
  with open(input_file, "r") as f:
    data = f.read()
  sets = data.split("[Results]")[1:]
  rets = {}
  for s in sets:
    s = s.split("[/Results]")[0]
    testsize = int(s.split(":")[0].strip())
    testres_us = float(s.split(":")[1].split("ns")[0].strip()) / 1000 #ns -> us
    if testtype not in rets:
      rets[testtype] = {}
    rets[testtype][testsize] = testres_us
  return rets
  
rets = odsl_parse_file(args.testteefile)
rets = rets | odsl_parse_sgx_files(args.testsgxfile)

names = [
   'Encryption'
  ,'Decryption'
  ,'Mov'
  # ,'CMov'
  ,'Swap'      # rename to 'DiskSwap'
  ,'OCallSwap' # rename to 'EnclaveSwap'
]
xs = []
ys = []
for name in names:
  # xs += [f"{name}-32"]
  # ys += [rets[name][32]/rets['Mov'][512]]
  # xs += [f"{name}-512"]
  # ys += [rets[name][512]/rets['Mov'][512]]
  renamedName = name
  if name == "Swap":
    renamedName = 'DiskSwap'
  if name == "OCallSwap":
    renamedName = 'OCall'
  xs += [f"{renamedName}-4096"]
  ys += [rets[name][4096]/rets['Mov'][4096]]

plt.barh(xs,ys,height=0.9,align="center")
axes = plt.gca()
axes.set_aspect(0.8)
plt.savefig(args.outfile, bbox_inches='tight')