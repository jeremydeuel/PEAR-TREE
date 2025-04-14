f = "efficiency_stats.txt"
# generate this file by running seff JOBID >> efficiency_stats.txt for all completed jobids
import pandas as pd
stats = {}

if __name__ == '__main__':
    current_job = None
    current_stats = {}
    with open(f, 'r') as f:
        for l in f:
            l = l.strip()
            if len(l):
                if l[:8] == 'Job ID: ':
                    if current_job is not None:
                        stats[current_job] = current_stats
                    current_job = l[8:]
                    current_stats = {}
                    continue
                if l[:7] == 'State: ':
                    current_stats['state'] = l[7:]
                if l[:14] == 'CPU Utilized: ':
                    current_stats['cpu'] = l[14:]
                if l[:17] == 'Memory Utilized: ':
                    current_stats['mem'] = l[17:]
    stats[current_job] = current_stats
    d = pd.DataFrame(stats).transpose()
    print(d)
    print(d.mem.max())
