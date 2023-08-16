#!/usr/bin/env python
import sys
import csv
import math
csv_file = sys.argv[1]
count = 10000000
threshold = 0
with open(csv_file,'r') as histo:
    csvreader = csv.reader(histo,delimiter=' ')
    for row in csvreader:
        if count >= int(row[-1]):
            count = int(row[-1])
            threshold = int(int(row[0])/2)
        else: #found local min
            if threshold<2:
                sys.exit(0)
            else:
                sys.stdout.write(str(threshold))
                sys.exit(0)
