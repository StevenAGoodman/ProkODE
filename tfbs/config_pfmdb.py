from Bio import motifs
import random

file_paths = ['regtransbase.meme','collectf.meme','dpinteract.meme','fan2020.meme','SwissRegulon_e_coli.meme']

mots = []
for file in file_paths:
    record = motifs.parse(open(f'tf_pfms/{file}', 'r'), 'minimal')
    for m in record:
        m.name = m.name[:m.name.find('_')]
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
