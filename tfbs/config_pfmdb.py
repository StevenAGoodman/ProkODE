from Bio import motifs
import random
import re
import os

directory = os.fsencode('tf_pfms_prok')
mots = []

for file in os.listdir(directory):
    filename = os.fsdecode(file)

    record = motifs.parse(open(f'tf_pfms_prok/{filename}', 'r'), 'minimal')
    for m in record:
        foo = re.search("(^| |_)([a-zA-Z][a-zA-Z][a-zA-Z][a-zA-Z]?)( |_|$)", m.name)
        m.name = foo.group().replace('_', '') if foo != None else 'abort'
        if m.name == 'abort':
            continue
        mots.append(m)

out = motifs.write(mots,'JASPAR')
open('pfmdb.txt', 'w').write(out)

## Prediction and Analysis of Transcription Factor Binding Sites: Practical Examples and Case Studies Using R Programming
# uniprobe

ms = []
record = motifs.parse(open('pfmdb.txt', 'r'), 'jaspar')
for m in record:
    m.matrix_id = f'MA{random.randint(0,9)}{random.randint(0,9)}{random.randint(0,9)}{random.randint(0,9)}.{random.randint(0,9)}'
    ms.append(m)

out = motifs.write(ms, 'JASPAR')
open('pfmdb.txt', 'w').write(out)
