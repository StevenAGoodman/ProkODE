import re
import pandas as pd

def annotation():
    with open("genomic.gbff","r") as file:
        out = "geneid\tstart\tend\tsynonyms\n"
        get_synonym = False
        spec = False
        prev_lin = 0
        for line in file:
            if spec:
                stri += " " + line.strip()[:-1]
                out += f"{stri.split("; ")}\n"
                spec = False
            if get_synonym:
                search2 = re.search('/gene_synonym=".+"', line)
                search3 = re.search('/gene_synonym="[^"]+\n', line)
                if search2 != None:
                    get_synonym = False
                    out += f"{search2.group()[15:-1].split("; ")}\n"
                elif search2 == None and search3 != None:
                    spec = True
                    stri = search3.group()[15:-1]


            search = re.search('/gene=".+"', line)
            if search != None:
                geneid = search.group()[7:-1]
                start = re.search("\s(complement\()?\d+\.\.\d+", prev_lin)
                if start == None:
                    continue
                start = start.group()
                end = start[start.find("..")+2:]
                start = start[12:start.find("..")] if start.find("nt(") >= 0 else start[1:start.find("..")]
                prev_lin = line
                out += f"{geneid}\t{start}\t{end}\t"
                get_synonym = True
                continue
            else:
                prev_lin = line
                continue


    out = out.split("\n")
    res = []
    foo_prev = "na"
    for lin in out:
        foo = lin[:lin.find("\t")]
        if foo == foo_prev:
            foo_prev = foo
            continue
        else:
            res.append(lin + "\n")
            foo_prev = foo
            continue

    open("../../src/inputs/annotation.tsv", "w").writelines(res)

def only_cds():
    out_loc = "cds_genomic.gbff"
    with open("genomic.gbff","r") as file:
        with open(out_loc, 'a') as outf:
            outf.write("""LOCUS       NC_000913            4641652 bp    DNA     circular CON 09-MAR-2022
DEFINITION  Escherichia coli str. K-12 substr. MG1655, complete genome.
ACCESSION   NC_000913
VERSION     NC_000913.3
DBLINK      BioProject: PRJNA57779
            BioSample: SAMN02604091
KEYWORDS    RefSeq.
SOURCE      Escherichia coli str. K-12 substr. MG1655
  ORGANISM  Escherichia coli str. K-12 substr. MG1655
            Bacteria; Pseudomonadota; Gammaproteobacteria; Enterobacterales;
            Enterobacteriaceae; Escherichia.
REFERENCE   1  (bases 1 to 4641652)
  AUTHORS   Riley,M., Abe,T., Arnaud,M.B., Berlyn,M.K., Blattner,F.R.,
            Chaudhuri,R.R., Glasner,J.D., Horiuchi,T., Keseler,I.M., Kosuge,T.,
            Mori,H., Perna,N.T., Plunkett,G. III, Rudd,K.E., Serres,M.H.,
            Thomas,G.H., Thomson,N.R., Wishart,D. and Wanner,B.L.
  TITLE     Escherichia coli K-12: a cooperatively developed annotation
            snapshot--2005
  JOURNAL   Nucleic Acids Res. 34 (1), 1-9 (2006)
   PUBMED   16397293
  REMARK    Publication Status: Online-Only
REFERENCE   2  (bases 1 to 4641652)
  AUTHORS   Hayashi,K., Morooka,N., Yamamoto,Y., Fujita,K., Isono,K., Choi,S.,
            Ohtsubo,E., Baba,T., Wanner,B.L., Mori,H. and Horiuchi,T.
  TITLE     Highly accurate genome sequences of Escherichia coli K-12 strains
            MG1655 and W3110
  JOURNAL   Mol. Syst. Biol. 2, 2006 (2006)
   PUBMED   16738553
REFERENCE   3  (bases 1 to 4641652)
  AUTHORS   Blattner,F.R., Plunkett,G. III, Bloch,C.A., Perna,N.T., Burland,V.,
            Riley,M., Collado-Vides,J., Glasner,J.D., Rode,C.K., Mayhew,G.F.,
            Gregor,J., Davis,N.W., Kirkpatrick,H.A., Goeden,M.A., Rose,D.J.,
            Mau,B. and Shao,Y.
  TITLE     The complete genome sequence of Escherichia coli K-12
  JOURNAL   Science 277 (5331), 1453-1462 (1997)
   PUBMED   9278503
REFERENCE   4  (bases 1 to 4641652)
  AUTHORS   Arnaud,M., Berlyn,M.K.B., Blattner,F.R., Galperin,M.Y.,
            Glasner,J.D., Horiuchi,T., Kosuge,T., Mori,H., Perna,N.T.,
            Plunkett,G. III, Riley,M., Rudd,K.E., Serres,M.H., Thomas,G.H. and
            Wanner,B.L.
  TITLE     Workshop on Annotation of Escherichia coli K-12
  JOURNAL   Unpublished
  REMARK    Woods Hole, Mass., on 14-18 November 2003 (sequence corrections)
REFERENCE   5  (bases 1 to 4641652)
  AUTHORS   Glasner,J.D., Perna,N.T., Plunkett,G. III, Anderson,B.D.,
            Bockhorst,J., Hu,J.C., Riley,M., Rudd,K.E. and Serres,M.H.
  TITLE     ASAP: Escherichia coli K-12 strain MG1655 version m56
  JOURNAL   Unpublished
  REMARK    ASAP download 10 June 2004 (annotation updates)
REFERENCE   6  (bases 1 to 4641652)
  AUTHORS   Hayashi,K., Morooka,N., Mori,H. and Horiuchi,T.
  TITLE     A more accurate sequence comparison between genomes of Escherichia
            coli K12 W3110 and MG1655 strains
  JOURNAL   Unpublished
  REMARK    GenBank accessions AG613214 to AG613378 (sequence corrections)
REFERENCE   7  (bases 1 to 4641652)
  AUTHORS   Perna,N.T.
  TITLE     Escherichia coli K-12 MG1655 yqiK-rfaE intergenic region, genomic
            sequence correction
  JOURNAL   Unpublished
  REMARK    GenBank accession AY605712 (sequence corrections)
REFERENCE   8  (bases 1 to 4641652)
  AUTHORS   Rudd,K.E.
  TITLE     A manual approach to accurate translation start site annotation: an
            E. coli K-12 case study
  JOURNAL   Unpublished
REFERENCE   9  (bases 1 to 4641652)
  CONSRTM   NCBI Genome Project
  TITLE     Direct Submission
  JOURNAL   Submitted (08-MAR-2022) National Center for Biotechnology
            Information, NIH, Bethesda, MD 20894, USA
REFERENCE   10 (bases 1 to 4641652)
  AUTHORS   Blattner,F.R. and Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (30-JUL-2014) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
  REMARK    Protein update by submitter
REFERENCE   11 (bases 1 to 4641652)
  AUTHORS   Blattner,F.R. and Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (15-NOV-2013) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
  REMARK    Protein update by submitter
REFERENCE   12 (bases 1 to 4641652)
  AUTHORS   Blattner,F.R. and Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (26-SEP-2013) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
  REMARK    Sequence update by submitter
REFERENCE   13 (bases 1 to 4641652)
  AUTHORS   Rudd,K.E.
  TITLE     Direct Submission
  JOURNAL   Submitted (06-FEB-2013) Department of Biochemistry and Molecular
            Biology, University of Miami Miller School of Medicine, 118 Gautier
            Bldg., Miami, FL 33136, USA
  REMARK    Sequence update by submitter
REFERENCE   14 (bases 1 to 4641652)
  AUTHORS   Rudd,K.E.
  TITLE     Direct Submission
  JOURNAL   Submitted (24-APR-2007) Department of Biochemistry and Molecular
            Biology, University of Miami Miller School of Medicine, 118 Gautier
            Bldg., Miami, FL 33136, USA
  REMARK    Annotation update from ecogene.org as a multi-database
            collaboration
REFERENCE   15 (bases 1 to 4641652)
  AUTHORS   Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (07-FEB-2006) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
  REMARK    Protein updates by submitter
REFERENCE   16 (bases 1 to 4641652)
  AUTHORS   Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (10-JUN-2004) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
  REMARK    Sequence update by submitter
REFERENCE   17 (bases 1 to 4641652)
  AUTHORS   Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (13-OCT-1998) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
REFERENCE   18 (bases 1 to 4641652)
  AUTHORS   Blattner,F.R. and Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (02-SEP-1997) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
REFERENCE   19 (bases 1 to 4641652)
  AUTHORS   Blattner,F.R. and Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (16-JAN-1997) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
COMMENT     PROVISIONAL REFSEQ: This record has not yet been subject to final
            NCBI review. The reference sequence is identical to U00096.
            
            On Nov 3, 2013 this sequence version replaced NC_000913.2.
            Changes to proteins and annotation made on March 7, 2022.  Current
            U00096 annotation updates are derived from EcoCyc
            https://ecocyc.org/.  Suggestions for updates can be sent to
            biocyc-support@ai.sri.com. These updates are being generated from a
            collaboration  that includes EcoCyc, the University of Wisconsin,
            UniProtKB/Swiss-Prot, and the National Center for Biotechnology
            Information (NCBI).
            COMPLETENESS: full length.
FEATURES             Location/Qualifiers
""")
            in_cds = False
            in_gene = False
            transcript_found = False
            for line in file:
                search1 = re.search("\s\s\s\s\sCDS\s\s\s\s\s\s\s\s\s\s\s\s\s", line)
                search12 = re.search('\s\s\s\s\sgene\s\s\s\s\s\s\s\s\s\s\s\s', line)
                if search1 != None:
                    in_cds = True
                elif search12 != None:
                    in_gene = True
                
                if in_gene:
                    outf.write(line.replace('gene', 'CDS'))
                    if re.search('                     /db_xref="GeneID:', line) != None:
                        in_gene = False

                if in_cds:
                    if re.search("misc_feature    ", line) != None or re.search("mobile_element  ", line):
                        in_cds = False
                        continue

                    outf.write(line)
                    search2 = re.search('                     /translation=".+"', line)
                    search4 = re.search('                     /translation=".+', line)
                    if search2 != None:
                        transcript_found = True
                        in_cds = False
                        transcript_found = False
                    elif search4 != None:
                        transcript_found = True

                    if transcript_found:
                        search3 = re.search('                     [A-Z]+?"', line)
                        search5 = re.search('                     [A-Z]?"', line)
                        if search3 != None or search5 != None:
                            in_cds = False
                            transcript_found = False
                            continue
                else:
                    continue

def operons():
    annotation_df = pd.read_csv("../../src/inputs/annotation.tsv", delimiter="\t")
    out = "operonid\tstart\tend\tgeneids\n"
    started = False
    gene_found = False
    genes = []
    with open("operons.txt", "r") as operonsfile:
        for line in operonsfile:
            search1 = re.search(".+\t", line)
            if search1 == None:
                if started:
                    
                    out += f"{operonid}\t{start}\t{endcoord}\t{genes}\n"
                    genes = []
                operonid = re.search("\d+\n", line).group().replace("\n", "")
                print(operonid)
                started = True
                gene_found = True
            elif search1 != None and started:
                try:
                    coords = re.search("\t\d+\t\d+", line).group().replace("\t",",")[1:]
                    startcoord = coords[:coords.find(",")]
                    endcoord = coords[coords.find(",")+1:]
        
                    if gene_found:
                        start = startcoord
                        gene_found = False
                except:
                    None
                genes.append(search1.group()[1:search1.group()[1:].find("\t")+1].replace("\t",""))


    open("../../src/inputs/operons.tsv", "w").writelines(out)
 
    with open("../../src/inputs/operons.tsv", 'r') as r:
        out = ''
        for line in r:
            print(line)
            out += line.replace("'", '"')
        print(out)
    with open("../../src/inputs/operons.tsv", 'w') as f:
        f.write(out)
                


only_cds()