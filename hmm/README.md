# hmm_introgression
This is to extend the original himmix model to low-coverage and damaged genomic data. I rewrote scripts so it's covenient to adapt to my project for low-coverage data.

countfile can be modified from inputs from orignial version to get 5 columns: <br /> 
chr pos ancestral\_allele phase(0/1) num\_of\_derived\_alleles(>0)

generate maskfile(.bed) from coverage file, including 3 columns: <br />
chr start\_coord end\_coord

usage: <br />
**train and decode** <br />
main.py gt\_mode -count\_file *countfile* -mut\_file *mutfile* -mask\_file *maskfile* -out *prefix* -data\_type modern -phased(if phase is not all 0)<br />
**annotate** <br />
main.py anno -vcf *vcffile* -sample *3Neanerthals\_2Denisovans* -group1 *3Neanderthals* -group2 *2Denisovans* -anno *prefix* -called *calledfile* -count\_file *countfile* -map\_file *recmap* -map Shared\_Map(or other map) -phased(if phase is not all 0) -strict(refin called file, from first matching pos to last matching pos)


Different from orignial version: 1.train parameters using only bins with data;2.start from the first bin with data to last bin with data;3. treat all chromosomes as independent Markov chains(forward backward step, and decoding); additionally, the way I annotate fragments by matching them to published archaics; However, the results largely agree with original version

The script shall be used to reporduce results from our publication. For formal investigation, please cite and refer to Skov et al., 2018 PLOS Genetics
