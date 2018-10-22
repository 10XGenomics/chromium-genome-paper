import tenkit.hdf5
import kitten.utils
import pysam
import pandas as pd
import pyfasta
import numpy
import os

import resource
resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

#read in csv file with metadata and paths to the coverage.h5 files for the experiments being compared
#colnames in the input csv are: id1,id2,id1_libtype,id2_libtype,cell_line,sequencer,id1_gb,id2_gb,id1_template_mass,id2_template_mass,id1_length_sweep,id1_h5_path,id2_h5_path
data = pd.read_csv("input_csv")
bam_fn = "path_2_any_hg19_aligned_bam"
bam = pysam.Samfile(bam_fn)
fasta = pyfasta.Fasta("path_2_hg19-2.2.0_genome.fa")

def make_bed_df(gain,ref):
    in_region = False
    ranges = []
    
    for (index, depth) in enumerate(gain):
        if index==(len(gain)-1) and in_region:
            end_region=index+1
            ranges.append([ref,start_region,end_region]) 
        if not(in_region) and depth:
            in_region = True
            start_region = index
        if in_region and not(depth):
            in_region = False
            end_region = index
            ranges.append([ref,start_region,end_region]) 
    
    return(ranges)

# Args are  first cov.h5, second cov.h5, and BAM file.
def cov_diff(cov1, cov2, bam, fasta):
    pro_cov1 = 0
    pro_cov2 = 0
    
    cov1_df = []
    cov2_df = []
    
    for ref in bam.references:
        r1 = tenkit.hdf5.read_data_frame_indexed(cov1, [(ref, 0, 300000000)], query_cols=["mapq30_coverage_deduped"])
        r2 = tenkit.hdf5.read_data_frame_indexed(cov2, [(ref, 0, 300000000)], query_cols=["mapq30_coverage_deduped"])
        
        r1_th = 5*r1.mapq30_coverage_deduped.mean()
        r2_th = 5*r2.mapq30_coverage_deduped.mean()
        
        _seq = fasta[ref]
        seq = numpy.array(_seq, dtype="S1")
        valid_base = (seq[:] != "N")
        
        #Generate a boolean list for coverage gain
        gain1 = (r1.mapq30_coverage_deduped.between(5, r1_th, inclusive=False) & (r2.mapq30_coverage_deduped <= 5) & valid_base)
        gain2 = (r2.mapq30_coverage_deduped.between(5, r2_th, inclusive=False) & (r1.mapq30_coverage_deduped <= 5) & valid_base)
        
        #Generate the data frame for the bed files
        cov1_df.extend(make_bed_df(gain1.tolist(),ref))
        cov2_df.extend(make_bed_df(gain2.tolist(),ref))
        
        #Sum the bases that meet criteria
        pro_cov1 = pro_cov1 + gain1.sum()
        pro_cov2 = pro_cov2 + gain2.sum()

        # Log the thresholds used for each chromosome
        print ref,r1_th,r2_th
        
    cov1_df = pd.DataFrame(cov1_df,columns=["chrom","start","end"])
    cov2_df = pd.DataFrame(cov2_df,columns=["chrom","start","end"])
    return (pro_cov1, pro_cov2,pro_cov1-pro_cov2),cov1_df,cov2_df

records = []
for index, row in data.iterrows():
    r,cov1_df,cov2_df = cov_diff(row['id1_h5_path'],row['id2_h5_path'],bam,fasta)
    records.append(r)
    print "Finished processing %d and %d: %s" %(row['Var1'],row['Var2'],r)
    if row['id1_length_sweep']=="No":
        fn_prefix = os.path.join('BED',row['cell_line']+'_'+str(row['id1'])+'_'+str(row['id2']))
        cov1_df.to_csv(fn_prefix+"_chromium.bed", encoding='utf-8', index=False,sep='\t',header=False)
        cov2_df.to_csv(fn_prefix+"_truseq.bed", encoding='utf-8', index=False,sep='\t',header=False)

index = ["Cov_in_10x","Cov_in_TruSeq","Delta"]
df = pd.DataFrame(records,columns=index)
result = pd.concat([data,df], axis=1, join_axes=[data.index])
result.to_csv("cov_gain_results.csv",encoding='utf-8', index=False)
