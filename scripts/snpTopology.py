from snakemake.script import snakemake
import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import statistics

from Bio import SeqIO
from Bio import Seq

import deconstructSigs

from scipy.stats import zscore

ALL_A_THALIANA_CHROMOSOMES = ["1", "2", "3", "4", "5", "Pt", "Mt"]
NUMBERED_A_THALIANA_CHROMOSOMES = ["1", "2", "3", "4", "5"]
A_THALIANA_CHROMOSOME_LENGTHS = {"1":30427671,
                               "2":19698289,
                               "3":23459830,
                               "4":18585056,
                               "5":26975502,
                               "Pt":154478,
                               "Mt":366924}
#centromere loci taken from chatgpt need to find source for this informaiton
A_THALIANA_CENTROMERE_LOCATIONS = { "1":[14150000,15100000], 
                                    "2":[3620000,4490000],
                                    "3":[13000000,14000000],
                                    "4":[1850000,2780000],
                                    "5":[10740000,11550000]}
SAMPLES = ['Sin1', 'Son1', 'San8', 'San11',]
NUCLEOTIDES = ["A","C","G","T"]
PURINES = ["A","G"]
PYRIMIDINES = ["C","T"]
COMPLIMENTS = str.maketrans("ATGC", "TACG")
SUBSTITUTIONS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

CHROMSOME_COLORS = {"1": '#0072BD',
                    "2": '#D95319',
                    "3": 'grey',
                    "4": '#EDB120',
                    "5": '#77AC30'}



def vcf_read_in(sampleID, in_vcf):
    """Reads .vcf file into chromosome separated dataframes, sotred in a dictionary
    (one entry in dict for each chromosome).

    Parameters:
    -----------
    sampleID : string
      the sample ID label on the vcf to be read in (vcf file should be stored in local directory)
    chromosome_list : list[string]
      the list of known chromosome names (formated to match the 'chr' field of vcf)

    Return:
    -------
    genome_dict : dictionary{string : pandas dataframe}
      a dictionary of chromosome variant dataframes

    Nested Functions:
    -----------------
    chrosomosome_separator(genome_df)
    extract_allele_frequency(genome_dict, sampleID)
    """
    global NUCLEOTIDES, ALL_A_THALIANA_CHROMOSOMES

    def chromosome_separator(genome_df):
        """Separates the variants in a whole genome variant dataframe into seperate data frames
        (as many dataframes as there chromosomes). Separated chromosome variant data frames are
        stored in a dictionary with key values for each chromosome.

        Parameters:
        -----------
        genome_df : pandas dataframe
            the variants called across the whole genome

        Return:
        -------
        genome_dict : dictionary{string : pandas dataframe}
            a dictionary of chromosome variant dataframes
        """
        global ALL_A_THALIANA_CHROMOSOMES
        nonlocal sampleID

        genome_dict = {}

        total = len(genome_df.index)
        assigned = 0

        for chromosome in ALL_A_THALIANA_CHROMOSOMES:
            genome_dict[chromosome] = genome_df[genome_df["chr"] == chromosome]
            genome_dict[chromosome] = genome_dict[chromosome].reset_index(drop=True)
            assigned += len(genome_dict[chromosome].index)
        
        print(f"{assigned}/{total} {sampleID} Variants Assigned")

        return genome_dict


    def extract_allele_frequency(genome_df):
        nonlocal sampleID

        AF_field = []

        for index, variant in genome_df.iterrows():
            format_field = variant["format"]
            sample_field = variant[sampleID]
            format_field_array = format_field.split(':')
            sample_field_array = sample_field.split(':')
            sample_format_dict = dict(zip(format_field_array, sample_field_array))
            AF_field.append(sample_format_dict['AF'])

        genome_df["AF"] = AF_field

        return genome_df

    whole_genome_df = pd.read_csv(in_vcf,
                        sep="\t",
                        comment="#",
                        dtype="str",
                        header=None,
                        names=["chr","pos","id","ref","alt","qual","filter","info","format", sampleID, "Ws2"])
    whole_genome_df["pos"] = pd.to_numeric(whole_genome_df["pos"])

    format_field = whole_genome_df["format"][0]
    format_field_array = format_field.split(':')
    if 'AF' in format_field_array:
      whole_genome_df = extract_allele_frequency(whole_genome_df)

    genome_dict = chromosome_separator(whole_genome_df)

    return genome_dict


class trinucSNP:

    def __init__(self, chromosome, position, ref, alt, fasta):

        def normalize_sbs(chromosome, position):
            global COMPLIMENTS, NUCLEOTIDES
            nonlocal fasta

            for base in self.trinucleotide_context:
                if base not in NUCLEOTIDES:
                    self.valid = False
                    print(f'Invalid nucleotide in trinucleotide context: {chromosome}:{position} {self.trinucleotide_context}')
            if self.ref in 'GA':
                self.trinucleotide_context = self.trinucleotide_context.translate(COMPLIMENTS)
                self.ref = ref.translate(COMPLIMENTS)
                self.alt = alt.translate(COMPLIMENTS)
        self.chromosome = chromosome
        self.position = position
        self.ref = ref
        self.alt = alt
        self.trinucleotide_context = str(fasta[chromosome-1].seq[(position-2):position+1]).upper()
        self.valid = True

        normalize_sbs(chromosome, position)

    def ssSNP_format(self):
        return "%s>%s" % (self.ref, self.alt)

    def triSNP_format(self):
        return "%s[%s]%s" % (self.trinucleotide_context[0], self.ssSNP_format(), self.trinucleotide_context[2])


class MutationalCatalog:

    def __init__(self, 
                 sampleID,
                 genome_var_dict, 
                 ref_fasta,
                 chromosome = None,
                 start_position = None,
                 end_position = None,
                 whole_genome = True):
        global SUBSTITUTIONS, NUCLEOTIDES
        
        def snp_iterator(variant_df):
            for index, variant in variant_df.iterrows():
                alts = variant['alt'].split(',')
                ref = variant['ref']
                for i in range(len(alts)):
                    if  (len(ref) == 1) and (len(alts[i]) == 1):
                        snp = trinucSNP(int(variant['chr']), int(variant['pos']), ref, alts[i], self.fasta)
                        if snp.valid:
                            self.sbsMatrix.loc[snp.triSNP_format()] += 1

        def whole_genome_catalog_propagation(genome_var_dict):

            for chromosome in NUMBERED_A_THALIANA_CHROMOSOMES:
                chromosome_df = genome_var_dict[chromosome]
                snp_iterator(chromosome_df)
        
        def windowed_catalog_propagation(genome_var_dict, chromosome, start_position, end_position):

            chromosome_df = genome_var_dict[chromosome]
            window_df = chromosome_df[(chromosome_df["pos"] >= start_position) & (chromosome_df["pos"] < end_position)]

            snp_iterator(window_df)

        self.fasta = list(SeqIO.parse(ref_fasta, 'fasta'))


        #generating empty matrix
        substitution_types = []
        for substitution in SUBSTITUTIONS:
            for context5prime in NUCLEOTIDES:
                for context3prime in NUCLEOTIDES:
                    substitution_types.append("%s[%s]%s" % (context5prime, substitution, context3prime))

        self.sbsMatrix = pd.DataFrame(index=substitution_types)
        self.sbsMatrix.index.name = "Type"

        #pre-filling matrix with 0's
        empty_sample_field = []
        for type in self.sbsMatrix.index.tolist():
            empty_sample_field.append(0)
        self.sbsMatrix[sampleID] = empty_sample_field

        if whole_genome:
            whole_genome_catalog_propagation(genome_var_dict)
        else:
            windowed_catalog_propagation(genome_var_dict, chromosome, start_position, end_position)

def sig11Contribution(mutational_catalog, sampleID):
    #tranlate mutational catalog matrix object into dictionary for input into deconstructSigs
    mutational_catalog_dict = {}
    for SNPtype in mutational_catalog.sbsMatrix.index:
        mutational_catalog_dict[SNPtype] = mutational_catalog.sbsMatrix.loc[SNPtype][sampleID]

    #construct DS model and extract signature weights
    DS_model = deconstructSigs.DeconstructSigs( context_counts= mutational_catalog_dict, 
                                                maf= './dummy/file/path')
    signatureWeights = DS_model.which_signatures()

    return signatureWeights[10]

def snpCountMatrixGenerator(genome_dict, chromosome, start_pos, end_position, AF_modifier):
    """Calulates the number of each type of SNP over a specified region.
    SNP counts are stored in a 3-dimensional dictionary(chromosome x reference x alternate).

    Parameters:
    -----------
    genome_dict : dictionary{(chromosome : string) : (variants : pandas dataframe)}
      a dataframe containing information all variants within region of interest
    chromosome : string
      the chromosome containing the region of interest
    start_pos, end_pos : int, int
      the start and end position of the region of interest, inclusive start, exclusive end


    Return:
    -------
    SNP_counts : dictionary{
        (reference : string) : dictionary{
            (alternate : string) : (count : int)}}
      a 2-dimensional dictionary of SNP counts split by reference and alternate identity
    """
    global NUCLEOTIDES, PURINES, PYRIMIDINES

    snpCountMatrix = {}
    for ref in NUCLEOTIDES:
      snpCountMatrix[ref] = {}
      for alt in NUCLEOTIDES:
        snpCountMatrix[ref][alt] = 0

    chromosome_df = genome_dict[chromosome]
    window_df = chromosome_df[(chromosome_df["pos"] >= start_pos) & (chromosome_df["pos"] < end_position)]
    for index, variant in window_df.iterrows():
      ref = variant["ref"]
      alts = variant["alt"].split(",")
      if AF_modifier:
        AFs = variant["AF"].split(",")
      for i in range(len(alts)):
        if len(ref) == 1 and len(alts[i]) == 1 and (ref in NUCLEOTIDES) and (alts[i] in NUCLEOTIDES):
          if AF_modifier:
              snpCountMatrix[ref][alts[i]] += float(AFs[i])
          else:
            snpCountMatrix[ref][alts[i]] += 1

    return snpCountMatrix


def tstv_calc(snpCountMatrix):
    """Calculates the ts/tv for a given from a given 2-dimensional dictionary of SNP counts.

    Parameters:
    -----------
    SNP_dict : dictionary{
        (reference : string) : dictionary{
            (alternate : string) : (count : int)}}
        a 2-dimensional dictionary of SNP counts split by reference and alternate identity

    Return:
    -------
    TS/TV : float
        the transition to transversion ratio for the given SNP dictionary representing the variation over a chosen region
    """
    global NUCLEOTIDES, PURINES, PYRIMIDINES

    TS = 0
    TV = 0

    for reference in NUCLEOTIDES:
        for alternate in NUCLEOTIDES:
            if (reference in PURINES and alternate in PYRIMIDINES) or (reference in PYRIMIDINES and alternate in PURINES):
                TV += snpCountMatrix[reference][alternate]
            else:
                TS += snpCountMatrix[reference][alternate]

    #do I want this to be included? how to properly handle zero divisor problems?
    if TV == 0:
        return TS

    return TS/TV

def CtoT_calc(snpCountMatrix):
    global NUCLEOTIDES, PURINES, PYRIMIDINES

    CtoT = 0
    altSNP = 0

    for reference in NUCLEOTIDES:
        for alternate in NUCLEOTIDES:
            #C>T change is indistinguishable from a G>A change since we 
            # dont know which strand changed first
            if (reference == "C" and alternate == "T") or (reference == "G" and alternate == "A"):
                CtoT += snpCountMatrix[reference][alternate]
            else:
                altSNP += snpCountMatrix[reference][alternate]

    if altSNP == 0:
        return 0

    return CtoT/altSNP


def whole_genome_tstv_calc(genome_dict):
    """Calculates the ts/tv over an entire given genome vairiant list given as a dictionary of dataframes.

    Parameters:
    -----------
    genome_dict : dictionary{(chromosome : string) : (variants : pandas dataframe)}
      a dataframe containing information all variants over the whole genome
    chromosome_list : list[string]
      the list of known chromosome names (formated to match the 'chr' field of vcf)

    Return:
    -------
    genome_tstv : float
      the transversion : transition ratio over the entire given genome variant list
    """
    nucs = ["A","C","G","T"]
    genome_SNP_counts = {}
    for ref in nucs:
      genome_SNP_counts[ref] = {}
      for alt in nucs:
        genome_SNP_counts[ref][alt] = 0


    for chromosome in genome_dict:
      chromosome_df = genome_dict[chromosome]
      chromosome_SNP_counts = snpCountMatrixGenerator(genome_dict, chromosome, 0, chromosome_df.iloc[-1]["pos"], AF_modifier=False)
      for reference in nucs:
        for alternate in nucs:
          genome_SNP_counts[reference][alternate] += chromosome_SNP_counts[reference][alternate]

    genome_tstv = tstv_calc(genome_SNP_counts)

    return genome_tstv



def collapse_topology(positions, topology):
    total_positions = []
    total_topology = []

    for chromosome in topology:
      total_topology = total_topology + topology[chromosome]
      for position in positions[chromosome]:
        total_positions.append(f'{chromosome}:{position}')
      total_positions.pop()
    
    return total_positions, total_topology

def get_max_ratio(topo):
    max_tstv = 0
    for chromosome in topo:
        for tstv_value in topo[chromosome]:
            if tstv_value > max_tstv:
                max_tstv = tstv_value
    return max_tstv

def get_threshold(positions, topology):

    total_positions, total_topology = collapse_topology(positions, topology)

    topo_stdev = statistics.stdev(total_topology)
    topo_avg = statistics.mean(total_topology)

    return ((2 * topo_stdev) + topo_avg)

def topology_mapping(sampleID, 
                     genome_dict,
                     method,
                     reference,
                     topology_plot,
                     distribution_plot,
                     window_size,
                     subsampling_ratio=1,
                     average_smoothing = False, average_pool_size = 10,
                     AF_modifier = False):
    """Generates ts/tv topography over the entirity of a given genome split by chromosomes
    (genome given as a variant dataframe split over chromosomes, stored in a dictionary)
    ts/tv is calculated over a given window size with a given subsampling ratio (default is 1 i.e., no subsampling).
    Mapping can be ploted as linear (with options for averaging smoothing and averaging bin sizing) or stair plot.
    Additional option for allele frquency modification of ts/tv calculation.

    Parameters:
    -----------
    genome_dict : dictionary{(chromosome : string) : (variants : pandas dataframe)}
      a dataframe containing information all variants over the whole genome
    plot_choise : string ('linear' or 'stair')
      the choice of plot for ts/tv topography (averaging and subsampling not supported for stair plot)
    window_size : int
      the size of the window over which to calculate ts/tv
    subsampling_ratio : float
      the step size for sliding window as a fraction of the window size (default is 1 i.e., no subsampling) (only applicable to linear plotting)
    average_smoothing : bool
      the choice of adjacent ts/tv data point smoothing (only applicable to linear plotting)
    average_pool_size : int
      the number of adjacent data points to average over (default is 10) (only applicable to linear plotting)
    AF_modifier : bool
      the choice of allele frequency modification in which AF of alternate is used instead of binary variant call (default is False)
    chromosome_lengths : dictionary{(chromosome : string) : (length : int)}
      the length of each chromosome in the genome (default is a_thaliana_chromosome_sizes)

    Return:
    -------
    None
    """
    global NUMBERED_A_THALIANA_CHROMOSOMES, A_THALIANA_CHROMOSOME_LENGTHS, A_THALIANA_CENTROMERE_LOCATIONS


    genome_topo = {}
    start_positions = {}
    for chromosome in NUMBERED_A_THALIANA_CHROMOSOMES:
    
        genome_topo[chromosome] = []
        start_positions[chromosome] = []

        for start_pos in range(1, A_THALIANA_CHROMOSOME_LENGTHS[chromosome], int(window_size*subsampling_ratio)):

            if method in ['C>T','TS/TV']:
              windowSNPMatrix = snpCountMatrixGenerator(genome_dict, chromosome, start_pos, start_pos + window_size, AF_modifier)
              if method == 'C>T':
                genome_topo[chromosome].append(CtoT_calc(windowSNPMatrix))
              elif method == 'TS/TV':
                genome_topo[chromosome].append(tstv_calc(windowSNPMatrix))
            elif method == 'Sig11':
                windowMutationalCatalog = MutationalCatalog(sampleID=sampleID, 
                                                            genome_var_dict=genome_dict, 
                                                            ref_fasta=reference, 
                                                            chromosome=chromosome, 
                                                            start_position=start_pos, 
                                                            end_position=start_pos + window_size, 
                                                            whole_genome=False)
                genome_topo[chromosome].append(sig11Contribution(
                                                        mutational_catalog= windowMutationalCatalog,
                                                        sampleID= sampleID))
            start_positions[chromosome].append(start_pos)

    
    peak_threshold = get_threshold(start_positions, genome_topo)


    def topology_stair_plotting(positions, topo):
        global A_THALIANA_CHROMOSOME_LENGTHS, A_THALIANA_CENTROMERE_LOCATIONS, CHROMSOME_COLORS
        nonlocal method, window_size, average_smoothing, average_pool_size, AF_modifier, sampleID, peak_threshold

        fig, axs = plt.subplots(5, sharex=True, figsize=(9,10))
        fig.suptitle("%s\nBin Size: %d | AF Modifier: %s" % (sampleID, window_size, AF_modifier), fontsize=18, fontweight="bold")

        for chromosome in CHROMSOME_COLORS:
          positions[chromosome].append(positions[chromosome][-1] + window_size)

          subplot_index = int(chromosome) - 1

          axs[subplot_index].stairs(topo[chromosome], edges=positions[chromosome], color=CHROMSOME_COLORS[chromosome], fill=True)
          axs[subplot_index].set_title("Chromosome %s" % (chromosome), fontsize=12)

          axs[subplot_index].axvspan(A_THALIANA_CENTROMERE_LOCATIONS[chromosome][0], A_THALIANA_CENTROMERE_LOCATIONS[chromosome][1], color='grey', alpha=0.5)
          axs[subplot_index].axhline(y=peak_threshold, color='red', linestyle='--')

        for ax in axs.flat:
          ax.set(xlabel='Position', ylabel=method, ylim=(0, get_max_ratio(topo)))
          ax.label_outer()

        
        fig.savefig(topology_plot)
    

    def ratio_distribution(topo):
        global NUMBERED_A_THALIANA_CHROMOSOMES
        nonlocal method, sampleID, peak_threshold

        method_color_dict = { 'TS/TV' : 'gold', 
                              'C>T' : 'red',
                              'Sig11' : 'turquoise'}
        
        all_ratio_array = []

        for chromosome in NUMBERED_A_THALIANA_CHROMOSOMES:
            all_ratio_array = all_ratio_array + topo[chromosome]

        plt.figure(figsize=(10,5))
        plt.hist(x=all_ratio_array, color=method_color_dict[method], edgecolor= 'black', bins=100)
        plt.axvline(x= peak_threshold, color='red', linestyle="--")
        plt.suptitle(sampleID, fontsize=14, fontweight="bold")
        plt.title("Bin Size: %d | AF Modifier: %s" % (window_size, AF_modifier), fontsize=9)
        plt.xlabel(method)
        plt.ylabel('Windows')

        plt.savefig(distribution_plot)
        plt.close()

    ratio_distribution(genome_topo)
    topology_stair_plotting(start_positions, genome_topo)

    
    return start_positions, genome_topo

def region_file_generator(positions, topology, window_size, region_file):
  global NUMBERED_A_THALIANA_CHROMOSOMES, A_THALIANA_CHROMOSOME_LENGTHS

  total_positions, total_topology = collapse_topology(positions, topology)

  z_scores = zscore(total_topology)

  pass_positions = []
  for i in range(len(z_scores)):
    if z_scores[i] > 2:
      pass_positions.append(total_positions[i])
  
  with open(region_file, 'w') as out_f:
    for peak_position in pass_positions:
      chromosome = peak_position[0]
      start_position = int(peak_position[2:])
      end_position = start_position + window_size
      if start_position > window_size:
        start_position = start_position - window_size
      if end_position < A_THALIANA_CHROMOSOME_LENGTHS[chromosome]:
        end_position = end_position + window_size
      out_f.write(f'{chromosome}\t{start_position}\t{end_position}\n')
  






###MAIN###
whole_genome_dict = vcf_read_in(
                      in_vcf= snakemake.input['artifact_filtered_variants'],
                      sampleID= snakemake.wildcards['sampleID'])

passed_afModifier = snakemake.wildcards['afModifier']

if snakemake.wildcards['method'] == 'TSTV':
  passed_method = 'TS/TV'
elif snakemake.wildcards['method'] == 'CtoT':
  passed_method = 'C>T'
elif snakemake.wildcards['method'] == 'Sig11':
  passed_method = 'Sig11'
  passed_afModifier = False
else:
  raise("Please select a valid topography calulation method. Options: TSTV or CtoT")

positions, topology = topology_mapping( sampleID= snakemake.wildcards['sampleID'],
                                        genome_dict= whole_genome_dict,
                                        method= passed_method,
                                        reference= snakemake.input['reference_genome'],
                                        topology_plot= snakemake.output['topology_plot'],
                                        distribution_plot= snakemake.output['distribution_plot'],
                                        #how do I handle plot writing with new snakemake directive??????
                                        #plot_directory= f'{arguments.plotPath}/{sampleID}/{sampleID}{arguments.method}_w{arguments.windowSize}af{int(passed_afModifier)}',
                                        window_size= int(snakemake.wildcards['windowSize']), 
                                        AF_modifier= passed_afModifier)

region_file_generator(positions= positions,
                      topology= topology,
                      window_size= int(snakemake.wildcards['windowSize']),
                      region_file= snakemake.output['peak_regions'])
        
        
  