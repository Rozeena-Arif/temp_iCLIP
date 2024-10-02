import gffutils
import pybedtools
import sys

def process_gtf_files(human_gtf_path, sinv_gtf_path, output_sinv_gtf_path):
    # Import human GTF file
    human_gtf = pybedtools.BedTool(human_gtf_path)
    
    # Import SINV GTF file and update seqnames
    sinv_gtf = pybedtools.BedTool(sinv_gtf_path)
    sinv_gtf = sinv_gtf.each(lambda x: setattr(x, 'chrom', 'SINV') or x).saveas(output_sinv_gtf_path)
    
    # Combine human and SINV GTF files
    combined_gtf = human_gtf.cat(sinv_gtf, postmerge=False)
    combined_gtf.saveas(output_sinv_gtf_path)

# Example usage with Nextflow parameters
if __name__ == "__main__":
    human_gtf_path = sys.argv[1]
    sinv_gtf_path = sys.argv[2]
    output_sinv_gtf_path = sys.argv[3]
    #combined_gtf_path = sys.argv[4]

    # Call the function with provided paths
    process_gtf_files(human_gtf_path, sinv_gtf_path, output_sinv_gtf_path)

