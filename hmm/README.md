# hmm_introgression
See Skov et al., 2018 Plos Genetics

This is to extend the original himmix model to low-coverage and damaged genomic data. I rewrote scripts so it's covenient to adapt to my project for low-coverage data. 

Different from orignial version: 1.train parameters using only bins with data;2.start from the first bin with data to last bin with data;3. treat all chromosomes as independent Markov chains(forward backward step, and decoding); additionally, the way I annotate fragments by matching them to published archaics; However, the results largely agree with original version

For formal investigation, please cite and refer to Skov et al., 2018
