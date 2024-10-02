process CONCAT {
	tag "File saved to : ${baseDir}/assets/fasta/HS.SINV.fa"

	publishDir (
		path: "${baseDir}/assets/fasta",
		mode: 'copy',
		overwrite: 'true'
	)
  
	input:
		path fastaDir

	output:
		path 'HS.SINV.fa', emit: merge_fasta

	script:
		"""
		if [ ! -f ${baseDir}/assets/fasta/HS.SINV.fa ]; then
			cat ${fastaDir}/*.fa > HS.SINV.fa
		else
			echo "HS.SINV.fa already exists, creating symbolic link."
			ln -s ${baseDir}/assets/fasta/HS.SINV.fa ./
		fi
		"""
}

process GFF_TO_GTF {
	tag  "GFF to GTF, file saved to: ${baseDir}/assets/annotation/SINV_annotation.gtf"

	conda 'envs/agat.yml'

	publishDir (
		path: "${baseDir}/assets/annotation/",
		mode: 'copy',
		overwrite: 'true'
	)
  
	input:
		path sinvGFF

	output:
		path 'SINV_annotation.gtf', emit: sinv_gtf

	script:
		"""
		agat_convert_sp_gff2gtf.pl \
		--gff ${sinvGFF} \
		--gtf_version 3 \
		--output SINV_annotation.gtf
		"""
}

process REFORMAT {
	tag "Reformat GTF, file saved to: ${baseDir}/assets/annotation/reformated.gtf"

	conda 'envs/r-base.yml'

	publishDir (
	path: "${baseDir}/assets/annotation/",
	mode: 'copy',
	overwrite: 'true'
	)

	input:
		path GTFFile

	output:
		path 'reformated.gtf', emit: reformat_gtf

	script:
		"""
		Rscript ${baseDir}/bin/reformat_gtf.r ${GTFFile} reformated.gtf
		"""
}

process MERGE_GTF {
	tag "Merge host and sinv GTF, file saved to : ${baseDir}/assets/annotation/HS.SINV.gtf"

	conda 'envs/merge_gtf.yml'

	publishDir (
		path: "${baseDir}/assets/annotation/",
		mode: 'copy',
		overwrite: 'true'
	)

	input:
		path humanGTF
		path sinvGTF

	output:
		path 'HS.SINV.gtf', emit: merged_gtf

	script:
		"""
		cat ${humanGTF} ${sinvGTF} > HS.SINV.gtf
	"""
}

process ADD_MISSING_GENES {
	tag "Add missing gene and transcript names, file saved to: ${baseDir}/assets/annotation/HS.SINV_gn_tn.gtf"

	conda 'envs/pygtftk.yml'

	publishDir (
		path: "${baseDir}/assets/annotation/",
		mode: 'copy',
		overwrite: 'true'
	)

	input:
		path GTFFile

	output:
		path 'pre_gene_summary.txt', emit: pre_gene_summary
		path 'HS.SINV_gn.gtf', emit: gn_gtf
		path 'HS.SINV_gn_tn.gtf', emit: gn_tn_gtf
		path 'post_gene_summary.txt', emit: post_gene_summary

	script:
		"""
		gtftk count_key_values -u -i ${GTFFile} -k gene_id,gene_name,transcript_id,transcript_name -o pre_gene_summary.txt
		gtftk merge_attr -i ${GTFFile} -k gene_name,gene_id -d gene_name -o HS.SINV_gn.gtf
		gtftk merge_attr -i HS.SINV_gn.gtf -k transcript_name,transcript_id -d transcript_name -o HS.SINV_gn_tn.gtf
		gtftk count_key_values -u -i ${baseDir}/assets/annotation/HS.SINV_gn.gtf -k gene_id,gene_name,transcript_id,transcript_name -o 'post_gene_summary.txt'
		"""
}

process FLATTEN_ANNOTATION {
	tag "Flatten the annotation, file saved to : ${baseDir}/assets/annotation/HS.SINV.attribute-fix.flattened.w50s20.maptoid.txt.gz"

	conda 'envs/htseq.yml'

	publishDir (
		path: "${baseDir}/assets/annotation/",
		mode: 'copy',
		overwrite: 'true'
	)

	input:
		path GTFFile
		path fastaFile

	output:
		path 'HS.SINV.flattened.annotation.txt.gz', emit: flatten_gz
		path 'HS.SINV.attribute-fix.flattened.w50s20.txt.gz', emit: sliding_window
		path 'HS.SINV.attribute-fix.flattened.w50s20.maptoid.txt.gz', emit: mapped_to_id
		path 'genome_size.txt', emit: genome_size

	script:
		"""
		htseq-clip annotation \
		-g ${GTFFile} \
		-u gene_id \
		-n gene_name \
		-t gene_biotype \
		--splitExons \
		--unsorted \
		-o HS.SINV.flattened.annotation.txt.gz

		htseq-clip createSlidingWindows \
		-i HS.SINV.flattened.annotation.txt.gz \
		-w 50 \
		-s 20 \
		-o HS.SINV.attribute-fix.flattened.w50s20.txt.gz

		htseq-clip mapToId \
		-a HS.SINV.attribute-fix.flattened.w50s20.txt.gz \
		-o HS.SINV.attribute-fix.flattened.w50s20.maptoid.txt.gz

		cut -f1,2 ${baseDir}/assets/fasta/HS.SINV.fa.fai > genome_size.txt
		"""
}

process REMOVE_PHIX_READS {
	tag "Removing PhiX reads, file saved to : ${params.rawFileDir}/PhiX/${sampleName}_rmPhiX.${params.fqExt}"

	conda 'envs/bbmap.yml'

	publishDir (
		path: "${params.rawFileDir}/PhiX/",
		mode: 'copy'
	)

	input:
		tuple val(sampleName)

	output:
		path "${sampleName}_rmPhiX.${params.fqExt}", emit: phix

	script:
		"""
		bbduk.sh \
		in="${params.rawFileDir}/${sampleName}.${params.fqExt}" \
		out="${sampleName}_rmPhiX.${params.fqExt}" \
		ref="${params.phixRef}"
		"""
}

process FASTQC {

	tag "FastQC of file : ${params.plots}/fastQC/${phixFastq}"

	conda 'envs/fastqc.yml'

	publishDir (
		path: "${params.plots}",
		mode: 'copy',
		overwrite: 'true'
		)

	input:
		path phixFastq

	output:
		path "fastQC/${phixFastq.baseName}", emit: fqc

	script:
		"""
		mkdir -p fastQC/${phixFastq.baseName}
		fastqc ${phixFastq} -o "fastQC/${phixFastq.baseName}"
		"""
}

process MULTIQC {

	tag "MultiQC, file saved to : ${params.plots}/multiQC"

	publishDir "${params.plots}/multiQC", mode: 'copy'

	conda 'envs/multiqc.yml'

	input:
		path fastqcFiles

	script:
		"""
		multiqc -f ${fastqcFiles} -o "${params.plots}/multiQC"
		"""
}

process COUNT_BARCODE_FREQUENCIES {

	tag "File saved to : ${params.plots}/barcoded/${phixFastq.baseName}_barcode_freqs.txt"

	publishDir (
	path: "${params.plots}/barcoded",
	mode: 'copy',
	overwrite: 'true'
	)

	input:
	path phixFastq

	output:
	path "${phixFastq.baseName}_barcode_freqs.txt"

	script:
	"""
	echo "Count barcode frequencies"

	zcat ${phixFastq} | awk -v umi1_len=5 -v exp_bc_len=6 '{if (FNR%4==2) print substr(\$1,(umi1_len+1),exp_bc_len)}' | sort | uniq -c | sort -k1,1rn > "${phixFastq.baseName}_barcode_freqs.txt"

	echo "Demultiplex reads"
	"""
}

process DEMULTIPLEXED {
	publishDir "Trimming/Demultiplexed", mode: 'copy'

  input:
  val phixFastq

  output:
  path "${phixFastq.baseName}_unassigned_jemultiplexer.fastq.gz", emit: fastq

  script:
  """
  java -jar /software/je_1.2/je_1.2_bundle.jar demultiplex \
  F1="${phixFastq}" \
  BF="${params.baseDir}/NOP2_barcodes_file.txt" \
  RCHAR=':' \
  O="./" \
  UF1="${phixFastq.baseName}_unassigned_jemultiplexer.fastq.gz" \
  M="${phixFastq.baseName}_merged_jemultiplexer_out_stats.txt" \
  FASTQ_FILE_EXTENSION=fastq
  """
}

process FASTQC_DEMULTIPLEXED {

	conda 'envs/fastqc.yml'

	publishDir (
		path: "${params.plots}/Demultiplexed",
		mode: 'copy',
		overwrite: 'true'
	)

	input:
		path demuxFq

	output:
		path "fastQC/${demuxFq.baseName}", emit: fqc

	script:
		"""
		mkdir -p fastQC/${demuxFq.baseName}
		fastqc ${demuxFq} -o "fastQC/${demuxFq.baseName}"
		"""
}

process MULTIQC_DEMULTIPLEXED {

	publishDir "${params.plots}/Demultiplexed/multiQC", mode: 'copy'
	
	conda 'envs/multiqc.yml'

	input:
	path fastqcFiles

	script:
	"""
	multiqc -f ${fastqcFiles} -o "${params.plots}/Demultiplexed/multiQC"
	"""
}

process CUTADAPT {

	tag "File saved to :  ${params.trimming}/Adapter/${sampleName}_trimmed.${params.fqExt}"

	publishDir "${params.trimming}/Adapter/", mode: 'copy'

	conda 'envs/multiqc.yml'

	input:
		path unassignedFq
		tuple val(sampleName), val(val)

	output:
		path "${sampleName}_trimmed.${params.fqExt}", emit:trimmedFq
		path "${sampleName}_trimmed.fastq.tooshort.gz", emit:tooShortFq 
		path "${sampleName}_trimmed.fastq", emit: unzipFq

	script:
		"""
		cutadapt "${unassignedFq}" \
		-a AGATCGGAAGAGCGGTTCAG \
		-j 4 -e 0.1 -O 1 --nextseq-trim 10 --minimum-length 15 \
		-o "${sampleName}_trimmed.${params.fqExt}" \
		--too-short-output "${sampleName}_trimmed.fastq.tooshort.gz" \
		
		gunzip -c "${sampleName}_trimmed.${params.fqExt}" > "${sampleName}_trimmed.fastq"

		"""
}

//check if it creates result or not
process MULTIQC_ADAPTER_TRIMMED {

	tag "File saved to : ${params.plots}/Adapter/multiQC"

	publishDir "${params.plots}/Adapter/multiQC", mode: 'copy'

	conda 'envs/multiqc.yml'

	input:
		path trimmedFq

	script:
	"""
	mkdir -p "${params.plots}/Adapter/multiQC"
	multiqc -f ${trimmedFq} -o "${params.plots}/Adapter/multiQC"
	"""
}

process ACCESS_RNA_CONTAMINATION {
	
	tag "File saved to : ${params.trimming}/rRNA/${sampleName}.minus.rRNA.${params.fqExt}"

	conda 'envs/bbmap.yml'

	publishDir (
		path: "${params.trimming}/rRNA/",
		mode: 'copy'
	)

	//publishDir (
	//	path: "${params.trimming}/Adapter/",
	//	mode: 'copy'
	//)


	input:
		path inputFq
		tuple val(sampleName), val(val)

	output:
		path "${sampleName}.minus.rRNA.${params.fqExt}", emit: fastq
		path "${sampleName}.rRNA.stat.txt", emit: stat_txt
		//path "${sampleName}_trimmed.fastq", emit: unzipFq
		
	script:
		"""
		bbduk.sh \
		-Xmx24G \
		in="${params.trimming}/Adapter/${sampleName}_trimmed.${params.fqExt}" \
		out="${sampleName}.minus.rRNA.${params.fqExt}" \
		ref="${params.Hsap}" \
		k=25 \
		stats="${sampleName}.rRNA.stat.txt" \
		hdist=1

		#gunzip -f "${params.trimming}/Adapter/${sampleName}_trimmed.${params.fqExt}" > "${sampleName}_trimmed.fastq"
		"""
}

process STAR_GENOME_INDEX {

	tag "Generating STAR genome index"

	//publishDir (
	//	path: "${params.fastaDir}",
	//	mode: 'copy'
  //)
	storeDir "${params.fastaDir}"

	input:
		path genomeFa   
		path gtfFile             

	output:
		path "STAR/SAindex", emit: genomeIndex
		path "STAR", emit: genomeDir
	
	script:
		"""
		STAR --runThreadN 2 --runMode genomeGenerate \
		--genomeFastaFiles "${genomeFa}" \
		--sjdbGTFfile "${gtfFile}" \
		--genomeDir "STAR" \
		--alignSJoverhangMin 8 \
		--sjdbOverhang 100 \
		--limitGenomeGenerateRAM 31000000000 \
		--genomeSAindexNbases 10 \
		--genomeSAsparseD 4
		"""
}

process STAR_GENOME_ALIGNMENT {

	tag "Aligning the reads to reference genome"

	storeDir "${params.results}/Alignment"

	//publishDir (
	//	path: "${params.results}/Alignment",
	//	mode: 'copy'
	//)

	input:
		path genomeDir
		path trimmedFq
		path gtfFile
		tuple val(sampleName)

	output:
		path "${sampleName}/Aligned.sortedByCoord.out.bam", emit: alignedFile

	script:
		"""
		STAR --runMode alignReads --runThreadN 2 \
		--outSAMtype BAM SortedByCoordinate \
		--genomeDir "${genomeDir}" \
		--sjdbGTFfile "${gtfFile}" \
		--readFilesIn "${trimmedFq}" \
		--outReadsUnmapped Fastx \
		--outFilterMismatchNmax 999 \
		--outFilterMultimapNmax 1 \
		--outFilterMismatchNoverLmax 0.04 \
		--outSJfilterReads Unique \
		--alignEndsType EndToEnd \
		--outFileNamePrefix "${sampleName}/"
		"""
}

process INDEX_ALIGNMENT {

	conda 'envs/samtools.yml'

	tag "Align the reads to reference genome"

	publishDir (
		path: "${params.results}/Alignment/",
		mode: 'copy'
	)

	input:
		path sortedAlignmentFile
		tuple val(sampleName)

	output:
		path "${sampleName}/Aligned.sortedByCoord.out.bam.bai", emit: bamIdx
		path "${sampleName}/idxstats.txt", emit: idxStat

	script:
	"""
	mkdir -p "${sampleName}"

	samtools index "${sortedAlignmentFile}" "${sampleName}/Aligned.sortedByCoord.out.bam.bai"

	samtools idxstats "${sortedAlignmentFile}" > "${sampleName}/idxstats.txt"
	"""
}

process MULTIQC_ALIGNED_BAM {

	tag "File saved to : ${params.results}/Alignment/multiQC"

	publishDir "${params.results}/Alignment/multiQC", mode: 'copy'

	conda 'envs/multiqc.yml'

	input:
		path idxStat

	script:
		"""
		multiqc -f ${idxStat} -o "${params.results}/Alignment/multiQC"
		"""
}

//format from here
process MARKDUPLICATE_ALIGNED_BAM {

	tag "File saved to : ${params.results}/Dedup/${sampleName}/${sampleName}_dedup.je.bam"
	
	publishDir "${params.results}/Dedup", mode: 'copy'

	input:
		path alignedBam
		tuple val(sampleName)

	output:
		path "${sampleName}/${sampleName}_dedup.je.bam", emit: dedupBam
		path "${sampleName}/${sampleName}_dedup.je.log", emit: dedupLog
		path "${sampleName}/${sampleName}_dedup.je.metrics", emit: dedupMetrics

	script:
	"""
	mkdir -p "${sampleName}"

	java -jar /software/je_1.2/je_1.2_bundle.jar markdupes \
	I="${params.results}/Alignment/${sampleName}/Aligned.sortedByCoord.out.bam" \
	O="${sampleName}/${sampleName}_dedup.je.bam" \
	M="${sampleName}/${sampleName}_dedup.je.metrics" \
	REMOVE_DUPLICATES=True \
	MM=1 > \
	"${sampleName}/${sampleName}_dedup.je.log"
	"""
}

process INDEX_DEDUP {

	tag "File saved to : ${params.results}/Dedup/${sampleName}/${sampleName}_dedup.je.bam.bai"

	conda 'envs/samtools.yml'

	tag "Samtools index the dedup BAM files"

	publishDir (
		path: "${params.results}/Dedup/",
		mode: 'copy'
	)

	input:
		path sortedAlignmentFile
		tuple val(sampleName)

	output:
		path "${sampleName}/${sampleName}_dedup.je.bam.bai", emit: bamIdx
		path "${sampleName}/${sampleName}_dedup.je.idxstats.txt", emit: idxStat

	script:
	"""
	mkdir -p "${sampleName}"

	samtools index "${sortedAlignmentFile}" "${sampleName}/${sampleName}_dedup.je.bam.bai"

	samtools idxstats "${sortedAlignmentFile}" > "${sampleName}/${sampleName}_dedup.je.idxstats.txt"
	"""
}

process MULTIQC_DEDUP_IDXSTATS {

	publishDir "${params.results}/Alignment/Dedup_stat_multiQC", mode: 'copy'

	conda 'envs/multiqc.yml'

	input:
		path idxStat

	script:
	"""
	multiqc -f ${idxStat} -o "${params.results}/Alignment/Dedup_stat_multiQC"
	"""
}


process RAW_READ_COUNT {

	tag "Read counts for trimmed, demultiplexed, aligned and dedup files"

	publishDir (
		path: "${params.trimming}/Adapter/",
		mode: 'copy'
	)

	publishDir (
		path: "${params.trimming}/Demultiplexed/",
		mode: 'copy'
	)

	publishDir (
		path: "${params.results}/Alignment/",
		mode: 'copy'
	)

	publishDir (
		path: "${params.results}/Dedup/",
		mode: 'copy'
	)

	input:
		path demultiplexedFq
		path adapterTrimmedFq
		path alignedBAM
		path dedupBAM

	output:
		path "${demultiplexedFq.baseName}.read_count.txt", emit: demuxCount
		path "${adapterTrimmedFq.baseName}.read_count.txt", emit: adapTrimCount
		path "${alignedBAM.baseName}.read_count.txt", emit: alignedCount
		path "${dedupBAM.baseName}.read_count.txt", emit: dedupCount

	script:
	"""
	count_dm=\$(zcat ${demultiplexedFq} | wc -l)
	reads_dm=\$((count_dm / 4))
	echo "${demultiplexedFq.baseName} \$reads_dm" > "${demultiplexedFq.baseName}.read_count.txt"

	count_at=\$(zcat ${adapterTrimmedFq} | wc -l)
	reads_at=\$((count_at / 4))
	echo "${adapterTrimmedFq.baseName} \$reads_at" > "${adapterTrimmedFq.baseName}.read_count.txt"
  
	count_alignedbam=\$(samtools view -c ${alignedBAM})
	echo "${alignedBAM.baseName} \$count_alignedbam" > "${alignedBAM.baseName}.read_count.txt"

	count_dedupbam=\$(samtools view -c ${dedupBAM})
	echo "${dedupBAM.baseName} \$count_dedupbam" > "${dedupBAM.baseName}.read_count.txt"
	"""
}

process XLSITE {

	tag "File saved to : ${params.results}/Xlsite/Collapsed/${sampleName}/${sampleName}.collapsed.bed"

	conda 'envs/bedtools.yml'

	tag "Bedtools shifted graph and collapsed file generation"

	publishDir (
		path: "${params.results}/Xlsite/Collapsed",
		mode: 'copy'
	)

	publishDir ( 
		path: "${params.results}/Xlsite/Shifted/bed",
		mode: 'copy'
	)

	input:
		path dedupBAM
		val genomeSize
		tuple val(sampleName)

	output:
		path "${sampleName}/${sampleName}.collapsed.bed", emit: collapsedBed
		path "${sampleName}/${sampleName}.shifted.bed", emit: shiftedBed

	script:
	"""
	mkdir -p "${sampleName}"

	bedtools bamtobed -i $dedupBAM > "${sampleName}/${sampleName}.collapsed.bed"

	bedtools shift -m 1 -p -1 -i "${sampleName}/${sampleName}.collapsed.bed" -g $genomeSize > "${sampleName}/${sampleName}.shifted.bed"
	"""
}

process BEDGRAPH_BED_FILE_GENERATION {

	conda 'envs/bedtools.yml'

	tag "Bedtools bedgraph"

	publishDir (
		path: "${params.results}/Xlsite/Shifted/bedgraph",
		mode: 'copy'
	)

	input:
		path shiftedBed
		val genomeSize
		val readCountFile

	output:
		path "${shiftedBed.baseName}.XL.RPM.+.bedgraph", emit: posBedGraph
		path "${shiftedBed.baseName}.XL.RPM.-.bedgraph", emit: negBedGraph

	script:
	"""
	base=\$(basename ${shiftedBed} .shifted.bed)
	TmpScale=\$(awk -v var="\$base" '\$1 == var { print \$2 }' ${readCountFile})
	scaleFactor=\$(bc <<< "scale=6; 1000000/\$TmpScale")

	outbgplus="${shiftedBed.baseName}.XL.RPM.+.bedgraph"
	outbgminus="${shiftedBed.baseName}.XL.RPM.-.bedgraph"

	bedtools genomecov -bg -strand + -5 -i ${shiftedBed} -g ${genomeSize} -scale \$scaleFactor > \$outbgplus
	bedtools genomecov -bg -strand - -5 -i ${shiftedBed} -g ${genomeSize} -scale \$scaleFactor > \$outbgminus
	"""
}

process BEDGRAPH_SORT {

	conda 'envs/bedtools.yml'

	tag "Bedtools bedgraph"

	publishDir (
		path: "${params.results}/Xlsite/Shifted/bedgraph",
		mode: 'copy'
	)

	input:
		path posBedGraph
		path negBedGraph

	output:
		path "${posBedGraph.baseName}.sorted.bedgraph", emit: sortedPosBedGraph
		path "${negBedGraph.baseName}.sorted.bedgraph", emit: sortedNegBedGraph

	script:
	"""
	sort -k1,1 -k2,2n ${posBedGraph} > "${posBedGraph.baseName}.sorted.bedgraph"
	sort -k1,1 -k2,2n ${negBedGraph} > "${negBedGraph.baseName}.sorted.bedgraph"
	"""
}

process BEDGRAPH_TO_BIGWIG {
	
	conda 'envs/bedgraphtobigwig.yml'

	tag "BedGraph to BigWig"

	publishDir (
		path: "${params.results}/Xlsite/Shifted/bw",
		mode: 'copy'
	)

	input:
		path posBedGraph
		path negBedGraph
		val genomeSize

	output:
		path "${posBedGraph.baseName}.bw", emit: bwPosBedGraph
		path "${negBedGraph.baseName}.bw", emit: bwNegBedGraph

	script:
	"""
	bedGraphToBigWig ${posBedGraph} ${genomeSize} "${posBedGraph.baseName}.bw"
	bedGraphToBigWig ${posBedGraph} ${genomeSize} "${negBedGraph.baseName}.bw"
	"""
}

