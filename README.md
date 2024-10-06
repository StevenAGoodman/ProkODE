# ProkODE: A Computational Model to Unpackage Time-Series Protein Concentrations during Prokaryotic Development

## Accessing Workspace

I am currently using external software. To allow others to edit this repo without the hassle of manually installing everything, I assume people have Visual Studio Code with the dev container functionality.

To work within this repo:

1. navigate to 'main' branch in this repository
2. download code to folder
3. open vscode
4. ensure 'Visual Studio Code Dev Containers' extension is installed
5. `ctrl + shift + p`: run 'Dev Containers: Open Folder in Container...'
6. select the unzipped folder you accessed through 'vscode_dev'

With this you should be able to edit the code with all the necessary packages. Please contact with any difficulties!

## Usage

<!--- in the future add a param for input folder --->

1. open in Docker (see above; right now it's just in vscode)
2. run `python ./src/run.py <prokode_directory> <reset_files_bool>`
   > - \<prokode_directory> <- base folder where ProkODE was intalled to
   > - \<reset_files_bool> <- if True, deletes all previously created files before running
3. run `python ./user_interface.py`

## Overview

Inputs:

- full genome sequence
- locations of annotated genes (indices within genome)
- args:
  - time frame to predict
  - list of proteins to show results

Outputs:

- protein concentrations across time frame (of listed proteins), visualizations and raw data

## Algorithm (deprecated)

1st layer = file; 2nd layer = function within file

> 1. /preprocessing/config_promo.py ->
>    > 1. writes to file (promoters.fa) containing sequences of the promoter regions of each gene
>    > 2. writes to file (promo_key.json) containing an ordered array of gene ids that corresponds to every 100 character interval in promoters.fa (i.e., the promoter of that gene)
> 2. /preprocessing/get_tfbs_decay.py
>    > 1. create_tfbs_file() -> writes file (../tfbs.csv) that contain rows for a transcription factor (tf), its target gene (tg), and its binding affinity (Kd) to tg
>    >    > 1. for gene in annotated genome:
>    >    >    > 1. search database of tf binding motifs (pscm_data.meme) for the tg id
>    >    >    >    > if not found, continue (as it is not a transcription factor)
>    >    >    > 2. scan all promoter regions in promoters.fa using external tool. for binding indices, divide by 100 (the len of each gene's promoter) and get coresponding gene id using promo_key.json
>    >    >    > 3. write to csv the tgs of each identified tf as well as the binding affinities, which are calculated in the external tool
>    > 1. create_decay_file() -> writes file (../decay_rates.csv) that contain rows for a gene and its mRNA decay rates (currently inputs random numbers; still searching for algorithm to predict decay rates...)
>    > 1. add_betas() -> adds a column to tfbs.csv for the beta value (the 'intrinsic regulatory parameter of the transcription factor'), currently random; still creating model that can predict values
> 3. /diffeq_plot.py
>    > 1. config_network_json() -> writes file (network.json) that can be processed into individual ODEs (in plot_system())
>    >    > 1. uses tfbs.csv, decay_rates.csv
>    >    > 2. see below for formating
>    > 2. plot_system() -> computes entire ode system and plots graphs of desired protein concentrations
>    >    > 1. using network.json, create differential equations for each gene based on effected genes (i will place the math in this readme soon)
>    >    > 2. plot graphs (haven't implemented the select proteins option yet)

/src/decay

1. filter motifs for only proteins that exist in genome
2. for each gene mrna/prot seq search for binding kd using ciiider (see if they're similar and can just be a constant)
3. return dict of [decaygene:[gene:kd,...]] pairs

## Testing Log

### Variables & testing notes for each
1. **binding probability functions**
   - P_i = N_i / (Nns * (N_i + Kd_i))
   - P_rnap = N_rnap / (Nns * (N_rnap + Kd_rnap-i)) for RNAP binding, \
   P_i = N_i / (N_i + Kd_i) for all else
   - SOME FUNCTION INCLUDING COLLISION PROBABILITY P_rnap = N_rnap / (Nns * (N_rnap + Kd_rnap-i)) for RNAP binding, \
   P_i = N_i / (N_i + Kd_i) for all else
2. **transcription function**
   - R_txn = B_all * P<sub>RNAP-gene-binding</sub> * R<sub>max_txn</sub>\
   s.t. R<sub>max_txn</sub> {mRNA / sec} = length_of_RNAP {nt} / transcription_rate {nt / sec}
2. **Beta function**
   - B_all = sum(P_i * B_i)
   - B_all = product(P_i * B_i)
      - if P_i = 0 for any i, it messes up the entire thing
   - B_all = sum(P_i * log10(B_i))
3. **translation function**
   - *Model A*: Translation_rate = P<sub>ribo-gene-binding</sub> * max_rate,\
   s.t. max_rate {prot / mRNA / sec} = length_of_ribo {aa} / translation_rate {aa / sec}
   - *Model B*: Translation_rate = P<sub>ribo-gene-binding</sub> * max_rate,\
   s.t. max_rate {prot / mRNA / sec} = length_of_ribo {aa} / translation_rate {aa / sec}
2. **mRNA decay model**
4. **protein decay model**
5. **RNAP change model**
6. **Ribosome change model**

### Testing & trial log

1. `1.1` *Configuration A*
   - variables:
      - model used:
   - trials:
      - trial 1
         - accuracy:
         - 
   - average accuracy:

2. `1.2`

### Overall validation log


## neat articles
https://iubmb.onlinelibrary.wiley.com/doi/epdf/10.1080/15216540400022441
https://en.wikipedia.org/wiki/Reaction_progress_kinetic_analysis
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5452729/
https://www.thoughtco.com/activation-energy-definition-ea-606348
https://chem.libretexts.org/Bookshelves/Introductory_Chemistry/Introductory_Chemistry/15%3A_Chemical_Equilibrium/15.02%3A_The_Rate_of_a_Chemical_Reaction
https://chem.libretexts.org/Bookshelves/Biological_Chemistry/Supplemental_Modules_(Biological_Chemistry)/Enzymes/Enzymatic_Kinetics/Michaelis-Menten_Kinetics

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6013389/#:~:text=Regulated%20proteolysis%20is%20a%20vital%20process%20that%20affects,specificity%20and%20selectivity%20in%20degrading%20substrates%20is%20key.
