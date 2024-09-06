# ProkODE: A Computational Model to Unpakage Time-Series Protein Concentrations during Prokaryotic Development

__NOTE__: I am currently using external software, so I just created a docker image for easier use...

## Overview

Inputs:
- full genome sequence
- locations of annotated genes (indices within genome)
- args:
  - time frame to predict
  - list of proteins to show results

Outputs:
- protein concentrations across time frame (of listed proteins), visualizations and raw data

## Algorithm

1st layer = file; 2nd layer = function within file

> 1. /preprocessing/config_promo.py ->
>> 1. writes to file (promoters.fa) containing sequences of the promoter regions of each gene
>> 2. writes to file (promo_key.json) containing an ordered array of gene ids that corresponds to every 100 character interval in promoters.fa (i.e., the promoter of that gene)
> 2. /preprocessing/get_tfbs_decay.py
>> 1. create_tfbs_file() -> writes file (../tfbs.csv) that contain rows for a transcription factor (tf), its target gene (tg), and its binding affinity (Kd) to tg
>>> 1. for gene in annotated genome:
>>>> 1. search database of tf binding motifs (pscm_data.meme) for the tg id
>>>>> if not found, continue (as it is not a transcription factor)
>>>> 2. scan all promoter regions in promoters.fa using external tool. for binding indices, divide by 100 (the len of each gene's promoter) and get coresponding gene id using promo_key.json
>>>> 3. write to csv the tgs of each identified tf as well as the binding affinities, which are calculated in the external tool
>> 1. create_decay_file() -> writes file (../decay_rates.csv) that contain rows for a gene and its mRNA decay rates (currently inputs random numbers; still searching for algorithm to predict decay rates...)
>> 2. add_betas() -> adds a column to tfbs.csv for the beta value (the 'intrinsic regulatory parameter of the transcription factor'), currently random; still creating model that can predict values
> 3. /diffeq_plot.py
>> 1. config_network_json() -> writes file (network.json) that can be processed into individual ODEs (in plot_system())
>>> 1. uses tfbs.csv, decay_rates.csv
>>> 2. see below for formating
>> 3. plot_system() -> computes entire ode system and plots graphs of desired protein concentrations
>>> 1. using network.json, create differential equations for each gene based on effected genes (i will place the math in this readme soon)
>>> 2. plot graphs (haven't implemented the select proteins option yet)

```{
    "tg1":[
        "tgdecay":float,
        "regulators":[
            "tf1": [
                "beta":float,
                "kd_tf":float,
            ]
            "tf2":[

            ]
            ...
        ]
    ],
    "tg2":[

    ],
    ...
} 
