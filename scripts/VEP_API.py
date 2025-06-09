from snakemake.script import snakemake

import requests, sys, json
import pandas as pd
from pprint import pprint
import time
import math

def vcf_read_in(sampleID, in_vcf):

    whole_genome_df = pd.read_csv(in_vcf,
                                  sep="\t", 
                                  comment="#", 
                                  dtype="str", 
                                  header=None, 
                                  names=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", sampleID, "Ws2"])

    return whole_genome_df

def vcf_write_out(tiered_df, out_vcf):
    with open(out_vcf, 'w') as out_f:
        out_f.write('#' + '\t'.join(tiered_df.columns.tolist()) + '\n')
        tiered_df.to_csv(out_vcf, sep="\t",index=False, header=False, mode='a')


def df_to_string_list(dataframe):
    return [' '.join(row) for row in dataframe[["CHROM","POS","ID","REF","ALT"]].values.tolist()]

        
def post_vcf_annotator(vcf_df):
    '''
    Annotates a vcf using the VEP post operation.
    Returns a df where there is one entry for each endpoint 
    (multiple entries in cases of variants that produced multiple endpoints).
    Currently does not include intergenic effects (VEP annotation of intergenic effects incompatible with other flagging and filtering)
    '''
    server = "https://rest.ensembl.org"
    request = "/vep/arabidopsis_thaliana/region"
    options = {"canonical" : 1}
    content_type='application/json'

    data_string_list = df_to_string_list(vcf_df)

    MAX_RETRIES = 8
    RETRY_DELAY = 5

    endpoints = []
    chunk_size = 50
    chunk = 1
    total_chunks = math.ceil(len(data_string_list)/chunk_size)

    for k in range(0, len(data_string_list), chunk_size): 
        data = {"variants" : data_string_list[k:k+chunk_size]}

        for attempt in range(MAX_RETRIES):
            try:
                my_r = requests.post(
                    server+request, 
                    headers= {  "Content-Type" :   content_type,
                                "Accept"       :   content_type },
                    json=data, 
                    params=options,
                    timeout=30)

                my_r.raise_for_status()

                endpoints = endpoints + my_r.json()
                break

            except requests.exceptions.RequestException as e:
                print(f'VCF Chunk: {chunk}/{total_chunks} Request Attempt {attempt+1} failed: {e}')
                if attempt < MAX_RETRIES - 1:
                    time.sleep(RETRY_DELAY)
                else:
                    raise
        chunk += 1
    
    annotated_df = pd.DataFrame(columns = vcf_df.columns.tolist() + ["VEP"])

    i = 0
    j = 0
    total = 0
    for endpoint in endpoints:
        if 'transcript_consequences' in endpoint:
            for effect in endpoint['transcript_consequences']:
                annotated_df.loc[j] = vcf_df.iloc[i].tolist() + [effect]
                j += 1
                total += 1
        # Option for including intergenic effects, not used for ease of use in other flags and filters
        # else:
        #     annotated_df.loc[j] = vcf_df.iloc[i].tolist() + [endpoint['intergenic_consequences']]
        #     j += 1
        i += 1

    return annotated_df

def canonical_flagging(annotated_df):
    '''
    Returns a dataframe including a flag field indicating whether the transcript affected is canonical.
    Values: bool (True | False)
    '''
    flagged_df = annotated_df.copy()
    canonical_field_list = []
    for index, variant in flagged_df.iterrows():
        if 'canonical' in variant['VEP']:
            canonical_field_list.append(True)
        else:
            canonical_field_list.append(False)

    flagged_df['CANONICAL'] = canonical_field_list

    return flagged_df

def consequence_flagging(annotated_df, 
                         coding_consequences = ['missense_variant', 
                                                'synonymous_variant', 
                                                'stop_gained', 
                                                'start_lost', 
                                                'splice_acceptor_variant']):
    '''
    Returns a dataframe including flag fields indicating transcript consequences and whether that consequence affects coding.
    Consequence Values: [str]
    Coding Values: bool (True | False)
    '''
    flagged_df = annotated_df.copy()
    consequence_field_list = []
    coding_field_list = []
    for index, variant in flagged_df.iterrows():
        consequences = []
        coding = False
        for consequence_term in variant['VEP']['consequence_terms']:
            consequences.append(consequence_term)
            if consequence_term in coding_consequences:
                coding = True
        consequence_field_list.append(consequences)
        coding_field_list.append(coding)
        
    flagged_df['CONSEQUENCE'] = consequence_field_list
    flagged_df['CODING'] = coding_field_list
    
    return flagged_df

def gene_panel_flagging(annotated_df, nicotinamide_gene_panel, circadian_gene_panel):
    '''
    Returns a dataframe including flag fields idciating the transcript affected and whether that transcript is included in a relevant gene panel.
    Transcript Values: str
    Gene Panel Values: bool (True | False)
    '''
    flagged_df = annotated_df.copy()

    nicotinamide_panel_df = pd.read_csv(nicotinamide_gene_panel, 
                                sep= "\t", 
                                dtype="str", 
                                header=0)
    nicotinamide_assoc_transcripts = set(nicotinamide_panel_df['Accession-1'])

    circadian_panel_df = pd.read_csv(circadian_gene_panel, 
                                sep= "\t", 
                                dtype="str", 
                                header=0)
    circadian_assoc_genes = set(circadian_panel_df['Gene stable ID'])    

    print(circadian_assoc_genes)
    
    gene_field_list = []
    transcript_field_list = []
    gene_panel_field_list = []
    for index, variant in flagged_df.iterrows():

        gene_panel_flags = []
        if variant['VEP']['transcript_id'] in nicotinamide_assoc_transcripts:
            gene_panel_flags.append('nicotinamide')
        if variant['VEP']['gene_id'] in circadian_assoc_genes:
            gene_panel_flags.append('circadian')
        if len(gene_panel_flags) == 0:
            gene_panel_flags.append('.')

        gene_panel_field_list.append(gene_panel_flags)

        gene_field_list.append(variant['VEP']['gene_id'])
        transcript_field_list.append(variant['VEP']['transcript_id'])
        
    flagged_df['GENE'] = gene_field_list
    flagged_df['TRANSCRIPT'] = transcript_field_list
    flagged_df['GENE_PANEL'] = gene_panel_field_list

    return flagged_df

def impact_flagging(annotated_df):
    flagged_df = annotated_df.copy()

    impact_field_list = []
    for index, variant in flagged_df.iterrows():
        impact_field_list.append(variant['VEP']['impact'])
    
    flagged_df['IMPACT'] = impact_field_list

    return flagged_df

def variant_tiering(annotated_df):
    tiered_df = annotated_df.copy()

    tier_field = []
    for index, variant in tiered_df.iterrows():
        if not variant['CODING'] or variant['IMPACT'] == 'MODIFER':
            tier_field.append(5)
        elif variant['GENE_PANEL'] != ['.']:
            tier_field.append(1)
        else:
            if variant['IMPACT'] == 'LOW':
                tier_field.append(4)
            elif variant['IMPACT'] == 'MODERATE':
                tier_field.append(3)
            elif variant['IMPACT'] == 'HIGH':
                tier_field.append(2)
    tiered_df['TIER'] = tier_field

    return tiered_df

# def gene_list_write_out(tiered_df, out_file):
#     with open(out_file, 'w') as o:
#         o.write('\n'.join(tiered_df['GENE']))




###MAIN###
sampleID = snakemake.wildcards['sampleID']

sampledf = vcf_read_in(
                sampleID= sampleID, 
                in_vcf= snakemake.input['region_filtered_variants']
                )

#Annotation
sampledf = post_vcf_annotator(sampledf)
print(f'{sampleID} After Annotation: {len(sampledf)}')



#Canonical Flagging and Filtering
sampledf = canonical_flagging(sampledf)
filteredsampledf = sampledf[sampledf['CANONICAL']].reset_index(drop=True)
filteredsampledf = filteredsampledf.drop(labels='CANONICAL', axis=1)
print(f'{sampleID} After Canonical Filter: {len(filteredsampledf)}')

#Flagging and Tiering
filteredsampledf = consequence_flagging(filteredsampledf)
filteredsampledf = impact_flagging(filteredsampledf)
filteredsampledf = gene_panel_flagging( annotated_df= filteredsampledf, 
                                        nicotinamide_gene_panel=snakemake.input['nicotinamide_gene_panel'],
                                        circadian_gene_panel=snakemake.input['circadian_gene_panel'])
filteredsampledf = filteredsampledf.drop(labels='VEP', axis=1)
tieredsampledf = variant_tiering(filteredsampledf)

#Tier Siphoning and Length Checks
tier5df = tieredsampledf[tieredsampledf['TIER'] == 5]
tier4df = tieredsampledf[tieredsampledf['TIER'] == 4]
tier3df = tieredsampledf[tieredsampledf['TIER'] == 3]
tier2df = tieredsampledf[tieredsampledf['TIER'] == 2]
tier1df = tieredsampledf[tieredsampledf['TIER'] == 1]
print(f'{sampleID} Tier 5: {len(tier5df)}')
print(f'{sampleID} Tier 4: {len(tier4df)}')
print(f'{sampleID} Tier 3: {len(tier3df)}')
print(f'{sampleID} Tier 2: {len(tier2df)}')
print(f'{sampleID} Tier 1: {len(tier1df)}')

print(tier1df.head())

#writing out list of involved genes (using official AT#G# format, not symbol)
# gene_list_write_out(tiered_df=tieredsampledf[tieredsampledf['TIER'] != 5], out_file= snakemake.output['gene_list'])

#writing out tiered variant calls as vcf (noncompliant with vcf standards)
vcf_write_out(
    tiered_df= tieredsampledf,
    out_vcf= snakemake.output['tiered_variants']
    )