file("_200seg7-10,103-106_log.txt",'w').writelines \
    ([' '.join([i.split()[ii] for ii in (3,5,7)]) + '\n' for i in file("200seg7-10,103-106_log.txt").readlines() \
      if i[0]=='['])
