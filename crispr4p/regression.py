#!/usr/local/bin/python

import datetime
import os
import time
starttime = time.time()

from crispr4p import PRECOMPUTED, AnnotationParser, PrimerDesign, \
    COORDINATES, SYNONIMS, FASTA



#names list
## from all synonims
names_list = [x[0] for x in AnnotationParser(COORDINATES, SYNONIMS).synonims_]

#run a full regression

mismatches_list = (0,)

os.system('rm -rf %s' % PRECOMPUTED)


#numer of mismatches
errors = []
pd = PrimerDesign(FASTA, COORDINATES, SYNONIMS, verbose=False, regression=True)
for ind, nmis in enumerate(mismatches_list):
    for inde, elem in enumerate(names_list):
        localArgs = ['--name', elem, '--mismatch', str(nmis)]
        _runtime = time.time()
        try:
            pd.runCL(localArgs)
        except:
            errors.append(' '.join(localArgs))
        print('Running regression:', inde, localArgs, time.time()-_runtime, len(errors))
if errors:
    with open('errors.txt', 'w') as fh:
        fh.write(str(datetime.datetime.now()) + '\n')
        fh.write('\n'.join(errors))

print('Running time:', time.time()-starttime)

