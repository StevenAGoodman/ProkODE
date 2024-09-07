from Bio import motifs

file_paths = ['regtransbase.meme','collectf.meme','dpinteract.meme','fan2020.meme','SwissRegulon_e_coli.meme']

mots = []
for file in file_paths:
    record = motifs.parse(open(f'tf_pfms/{file}', 'r'), 'minimal')
    for m in record:
        print(m.name)
        mots.append(m)

out = motifs.write(mots,'JASPAR')
open('pfmdb.jaspar', 'w').write(out)
