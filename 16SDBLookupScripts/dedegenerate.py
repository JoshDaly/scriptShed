#!/usr/bin/env python

def translateDegenerate(primer):
    ## this is crap I know, but does the job
    degenerates={"R": ["A", "G"],
                 "Y": ["C", "T"],
                 "S": ["G", "C"],
                 "W": ["A", "T"],
                 "K": ["G", "T"],
                 "M": ["A", "C"],
                 "B": ["C", "G", "T"],
                 "D": ["A", "G", "T"],
                 "H": ["A", "C", "T"],
                 "V": ["A", "C", "G"],
                 "N": ["A", "C", "G", "T"]}
    all=[]

    for nucl in primer:
        if nucl not in degenerates:
            if len(all) > 0:
                for idx,i in enumerate(all):
                    all[idx]=i+nucl
            else:
                all.append(nucl)
        else:
            p =  degenerates[nucl]
            if len(all)>0:
              all = all*len(p)
              x= len(all)/len(p)
              app_nucl= [p[i//x] for i in range(len(p)*x)]
              for idx, i in enumerate(all):
                  all[idx]=i+app_nucl[idx]
            else:
              all=p
    for p in set(all):
      yield p

import sys
for p in translateDegenerate(sys.argv[1]):
  print p
