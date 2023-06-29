import glob
import sys
filepath = sys.argv[1]
txt = glob.glob(filepath)
for filename in txt:
 with open(filename, 'r') as f:
  print('# ' + filename)
  lines = f.readlines()
  if len(lines) == 0:
    continue
  last = lines[-1]
  print('[', end='')
  for line in lines:
    time = line.split()[2]
    if line is last:
        print(time, end=']\n')
    else:
        print(time, end=', ')