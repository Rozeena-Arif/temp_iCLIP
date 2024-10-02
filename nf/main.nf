include { CONCAT } from './modules/genome_prep.nf'
include { GFF_TO_GTF } from './modules/genome_prep.nf'
include { REFORMAT } from './modules/genome_prep.nf'
include { MERGE_GTF } from './modules/genome_prep.nf'
include { ADD_MISSING_GENES } from './modules/genome_prep.nf'
include { FLATTEN_ANNOTATION } from './modules/genome_prep.nf'
include { REMOVE_PHIX_READS } from './modules/genome_prep.nf'
include { FASTQC } from './modules/genome_prep.nf'
include { MULTIQC } from './modules/genome_prep.nf'
include { COUNT_BARCODE_FREQUENCIES } from './modules/genome_prep.nf'
include { DEMULTIPLEXED } from './modules/genome_prep.nf'
include { FASTQC_DEMULTIPLEXED } from './modules/genome_prep.nf'
include { MULTIQC_DEMULTIPLEXED } from './modules/genome_prep.nf'
include { CUTADAPT } from './modules/genome_prep.nf'
include { MULTIQC_ADAPTER_TRIMMED } from './modules/genome_prep.nf'
include { ACCESS_RNA_CONTAMINATION } from './modules/genome_prep.nf'
include { STAR_GENOME_INDEX } from './modules/genome_prep.nf'
include { STAR_GENOME_ALIGNMENT } from './modules/genome_prep.nf'
include { INDEX_ALIGNMENT } from './modules/genome_prep.nf'
include { MULTIQC_ALIGNED_BAM } from './modules/genome_prep.nf'
include { MARKDUPLICATE_ALIGNED_BAM } from './modules/genome_prep.nf'
include { INDEX_DEDUP } from './modules/genome_prep.nf'
include { MULTIQC_DEDUP_IDXSTATS } from './modules/genome_prep.nf'
include { RAW_READ_COUNT } from './modules/genome_prep.nf'
include { XLSITE } from './modules/genome_prep.nf'
include { BEDGRAPH_BED_FILE_GENERATION } from './modules/genome_prep.nf'
include { BEDGRAPH_SORT } from './modules/genome_prep.nf'
include { BEDGRAPH_TO_BIGWIG } from './modules/genome_prep.nf'

def res_dir = new File("${params.plots}/fastQC")
if (!res_dir.exists()) {
        res_dir.mkdirs()
	}

metaFile = "${params.metaSheet}"

fq_channel = channel
  .fromPath(metaFile)
  .splitCsv(header: true, sep: ",") 
  .map { row -> tuple(row.sampleName, row.val) } 
  .view { it -> "Sample Name: ${it[0]}" }

println("Output file		:" + metaFile);

workflow {
  CONCAT(fastaDir=params.fastaDir)
  GFF_TO_GTF(sinvGFF=params.sinvGFF)
  REFORMAT(GTFFile=GFF_TO_GTF.out.sinv_gtf)
  MERGE_GTF(humanGTF=params.humanGTF, sinvGTF=REFORMAT.out.reformat_gtf)
  ADD_MISSING_GENES(GTFFile=MERGE_GTF.out.merged_gtf)
  FLATTEN_ANNOTATION(GTFFile=ADD_MISSING_GENES.out.gn_tn_gtf, fastaFile="${params.fastaDir}/HS.SINV.fa")
	REMOVE_PHIX_READS(sampleName=fq_channel)
	FASTQC(phixFastq=REMOVE_PHIX_READS.out.phix)
	MULTIQC(fastqcDir=FASTQC.out.fqc.collect())
	COUNT_BARCODE_FREQUENCIES(phixFasta=REMOVE_PHIX_READS.out.phix)
	DEMULTIPLEXED(phixFq=REMOVE_PHIX_READS.out.phix)
	FASTQC_DEMULTIPLEXED(demuxFq=DEMULTIPLEXED.out.fastq)
	MULTIQC_DEMULTIPLEXED(fastqcDir=FASTQC_DEMULTIPLEXED.out.fqc.collect())
	CUTADAPT(unassignedFq=DEMULTIPLEXED.out.fastq, fq_channel)
	////below is pending to check, output is not saving 
	MULTIQC_ADAPTER_TRIMMED(trimmedFq=CUTADAPT.out.trimmedFq.collect())

	ACCESS_RNA_CONTAMINATION(CUTADAPT.out.trimmedFq, fq_channel)
	STAR_GENOME_INDEX(genomeFa=CONCAT.out.merge_fasta, gtfFile=ADD_MISSING_GENES.out.gn_tn_gtf)
  STAR_GENOME_ALIGNMENT(genomeDir=STAR_GENOME_INDEX.out.genomeDir, trimmedFq=CUTADAPT.out.unzipFq, gtfFile=ADD_MISSING_GENES.out.gn_tn_gtf, fq_channel)
	INDEX_ALIGNMENT(sortedAlignmentFile=STAR_GENOME_ALIGNMENT.out.alignedFile, fq_channel)
	MULTIQC_ALIGNED_BAM(idxStat=INDEX_ALIGNMENT.out.idxStat)
	MARKDUPLICATE_ALIGNED_BAM(alignedBam=INDEX_ALIGNMENT.out.bamIdx, fq_channel)
	//INDEX_DEDUP(dedupBam=MARKDUPLICATE_ALIGNED_BAM.out.dedupBam, fq_channel)
	//MULTIQC_DEDUP_IDXSTATS(idxStat=INDEX_DEDUP.out.idxStat)
	//RAW_READ_COUNT(demultiplexedFq=DEMULTIPLEXED.out.fastq, adapterTrimmedFq=CUTADAPT.out.trimmedFq, alignedBAM=STAR_GENOME_ALIGNMENT.out.alignedFile, dedupBAM=MARKDUPLICATE_ALIGNED_BAM.out.dedupBam)
	//XLSITE(dedupBAM=MARKDUPLICATE_ALIGNED_BAM.out.dedupBam, genomeSize="${params.baseDir}/assets/fasta/sizes.genome", fq_channel)
	//BEDGRAPH_BED_FILE_GENERATION(shiftedBed=XLSITE.out.shiftedBed, genomeSize="${params.baseDir}/assets/fasta/sizes.genome", readCountFile=RAW_READ_COUNT.out.dedupCount)
	//BEDGRAPH_SORT(posBedGraph=BEDGRAPH_BED_FILE_GENERATION.out.posBedGraph,negBedGraph=BEDGRAPH_BED_FILE_GENERATION.out.negBedGraph)
	//BEDGRAPH_TO_BIGWIG(posBedGraph=BEDGRAPH_SORT.out.sortedPosBedGraph, negBedGraph=BEDGRAPH_SORT.out.sortedNegBedGraph, genomeSize="${params.baseDir}/assets/fasta/sizes.genome")
}

