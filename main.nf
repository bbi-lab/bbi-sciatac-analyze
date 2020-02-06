/* Default configuration parameters.
** Typically, many of these are given in a configuration file.
** The configuration file is specified as a parameter on the
** NextFlow command, for example,
**   nextflow run main.nf -c nextflow_parameters.config
**
** Run in grid qlogin shell with at least 16G of memory (my default allocation)
** Suggested command line parameters to use when running this script
**   nextflow run main.nf -w <work_dirname> -with-report <report_filename> -with-trace <trace_filename> -with-timeline <timeline_filename>
**
** Notes on NextFlow
**   o  no groovy/java code allowed in process{} block (apparently)
**      except in the script (multi-line) string where (at least)
**      groovy/java variables can be used. (It looks like println
**      statements may appear after the 'script:' statement and before
**      the first triple quotation marks. However, it appears that
**      this is an exception (it may be a NextFlow println rather
**      than a Groovy println function.
**   o  repeating myself, no groovy/java variables that are defined
**      in the 'input:' block are allowed (apparently) outside the
**      script string of a process{} block. However, global NextFlow
**      (and maybe Groovy) variables are allowed in at least some
**      non-script places in the process {} block, for example, in
**      directives.     
**   o  NextFlow 'file objects' are groovy/java 'path objects',
**      and they inherit the path object methods.
**   o  groovy code is allowed outside process{} blocks and other
**      special purpose NextFlow blocks (mostly, at least). This
**      may require some experimentation.
**   o  Hannah points out that NextFlow has a dump-hashes option, which
**      causes NextFlow to dump directory/file hashes. It uses these
**      hashes when it resumes so one can use this to help diagnose
**      resume problems.
**   o  NextFlow 'consumes' file type content in channels. This means
**      that a channel can be used only once in a script. In order to
**      use a channel's file content in more than one process block, the
**      channel must be duplicated using the .into{} operator. In
**      contrast to the file content, the 'val' content is not
**      consumed so it can be used in more than one process{} block.
**   o  publishDir: files in the pattern statement must be in an output channel
**   o  take care to escape dollar signs in shell and script (such as awk script)
**      substitutions. These unescaped dollar signs are interpreted as introducing
**      NextFlow/groovy variables in strings. Unfortunately, NextFlow mis-identifies
**      the problematic variable and its location in the code, which complicates
**      greatly finding the error. As an example, which confused me greatly, is
**
**          bedtools slop -i ${inTssRegions} -g ${inChromosomeSizes} -b {flanking_distance_per_base_tss_region}  \
**          | bedtools coverage -sorted -d -a stdin -b ${inTranspositionSites} \
**          | awk '{{ if (\$8 > 0) print \$0 }}' \
**          | gzip > \${TEMP_OUTPUT_FILE}
**
**      I ommitted the backslashes in the awk script's '$8' and '$0'. NextFlow
**      noted an error at 'stdin -b ${inTranspositionSites}', and it said that
**      the dollar sign must be escaped if it refers to a shell variable or the
**      'inTranspositionSites' must be enclosed in braces, which it is... I had
**      difficulty finding this error partly because the reported error location
**      occurs before the real error (odd in my experience), and the dollar signs
**      are technically not in the shell script, which seems a bit like a
**      rationalization.
**  o  NextFlow does not respect comments starting with '#' in the script block,
**     which may make sense. However, it complicates debugging because I want to
**     comment out portions of the script that I suspect may cause a NextFlow
**     error. (I have the impression that Groovy/Java comments are not allowed in
**     the script string also, which makes sense.)
**
** Notes on groovy
**   o  scope (these may not be entirely accurate):
**        o  any variable set without a type is a global variable, e.g., 'x = 3'. This
**           can be a source of problems if 'def' is omitted erroneously when defining
**           a list or dictionary in a function that is called more than once.
**        o  'def' or any type creates a local variable, e.g. 'int x' or 'def x'.
**        o  a local variable defined in the main script is accessible only within the main
**           script (not in functions defined within the main script, it appears!)
**        o  a local variable defined in a function is accessible only within that function
**        o  it appears that a variable 'def'ed in a code block, for example,
**           a loop or conditional block, is local to that block.
**        o  stack trace command: org.codehaus.groovy.runtime.StackTraceUtils.sanitize(new Exception()).printStackTrace()
**             or
**                def methodZ() {
**                def marker = new Throwable()
**                StackTraceUtils.sanitize(marker).stackTrace.eachWithIndex { e, i ->
**                   println "> $i ${e.toString().padRight(30)} ${e.methodName}"
**               }
**
** Notes on this script
**  o  I have tried to keep the comments and diagnostic output accurate
**     but I am certain that there are instances in which I cut and pasted
**     and forgot to edit the strings, or I edited the code and some
**     comments/strings no longer describe accurately the resulting code.
**  o  I regret the often long and complicated names; however,
**     I have not noticed a more satisfying alternative. I considered
**     abbreviations but strings of abbreviations in a name may
**     create more ambiguity and difficulty in comprehension. Part
**     of the challenge is distinguishing processes/variables/functions
**     and instances of these that are similar but not identical, of
**     which that are some in this script.
**  o  process block naming conventions
**        o  name ends with 'Process'
**  o  variable naming conventions
**       o  I have not done well establishing conventions
**          for variable names.
**  o  groovy function naming conventions
**       o  functions that set up a channel for input to a
**          process{} block end with 'ChannelSetup'; unless, the
**          process{} block has more than one input channel, in
**          which case the name ends with 'ChannelSetup<channel_content>'
**  o  NextFlow channel naming conventions
**       o  input channels to process block end with 'InChannel' unless
**          the process block has more than on input channel, in which
**          case the name ends with 'InChannel<channel_content>'
**       o  output channels from process block end with 'OutChannel' unless
**          the process block has more than on output channel, in which
**          case the name ends with 'OutChannel<channel_content>'
**       o  when a channel is copied (with an .into {} operator), I
**          append *Copy<nn> to the channel name, for example,
**          'getUniqueFragmentsOutChannelTranspositionSitesCopy01' (this is a
**          ungainly name, for sure).
**       o  in processes that have more than one input channel, I add
**          the a description of the file type to the end of the channel
**          name. I may append the file description to single file type
**          channels whose names do not adequately suggest the type of
**          file in it.
**       o  in processes that have more than one output channel, I add
**          the a description of the file type to the end of the channel
**          name. I may append the file description to single file type
**          channels whose names do not adequately suggest the type of
**          file in it.
**  o  use groovy/java functions to form InChannels from OutChannels
**     because
**       o  the groovy and java are well documented with good tutorials
**          and discussions
**       o  it may be easier to understand and follow groovy/java in
**          comparison to NextFlow operators for a few reasons such as
**          using intermediate variables and the ability to insert
**          diagnostic output statements. The function is generally
**          self-contained whereas the NextFlow operators are strung
**          together, sometimes in different locations inside and
**          outside the process blocks.
**       o  I believe that it is relatively easy to transcribe a
**          function into a stand-alone NextFlow or Groovy script for
**          testing.
**  o  in 'input' process{} statements, I try to keep each 'file' object as
**     a distinct (tuple) element that can be identified with the file() operator
**     within the statement. I use tuples for this purpose. The alternative is
**     to assign the file object as a map element but doing this may prevent one
**     from using the 'file' operator. It may not matter functionally. As an
**     aside, there are instances where I assign file names to map elements but
**     I do this only because the file names are strings rather than file
**     objects. (The file object is created by NextFlow for channel input and
**     output. The file name string is not (yet) linked to a channel.)
**
**  o  I am inclined to place a single file type in a process output channel.
**     Currently, there a few places in which there is a file and its index
**     file in a channel, and this may work out OK...I may later decide to place
**     the index files in a separate channel.
**  o  I dislike that there is more than one instance of the 'root' filename literal
**     string for most of the (input/output) files. Perhaps I can define groovy
**     variables to store these string literals: I am not certain that it can work.
**
** Tasks:
**   o  work on src/summarize_cell_calls.R and generate_cistopic_model.R to
**      remove dependencies on specific genome versions
**   o  add/update genomes
**   o  add dashboards/statistics
**   o  add plots for Cailyn
**   o  compact functions where possible
**   o  write program for aggregating runs
**   o  add in remaining downstream processes
**   o  consider using string variables for filenames
**   o  update file_map.docx
**   o  write JSON args.json file (perhaps write a sample-specific JSON file to each sample directory) (done)
**   o  fix reads_threshold in call_cells process (done)
**   o  store files in sub-directories (done)
**   o  build for new Shendure processors (no benefit)
**   o  clean repositories (done for now)
**
** Check:
**   o  check various global variable values such as flanking_distance*
**   o  check with all condition (parameters and genome files and ...) combinations
**   o  compare each process block to each easygrid pipeline block
**   o  compare outputs to those from Andrew's pipeline (done, for now)
**   o  run on BBI multi-sample run (done)
**   o  check for 'defs' where required (done, for now)
**   o  check booleans: e.g.,
**        isBarnyard: summarizeCellCallsSetupTestBarnyard
**        hasGeneScoreBed: hasGeneScoreBed
**   o  consider writing functions to replace often repeated code
**      such as checking for files in Groovy functions.
**   o  check diagnostic strings (e.g., asserts)
**   o  check comments for accuracy
**   o  check conditional code blocks
**   o  test run aggregation
**   o  check makePromoterSumIntervalsProcess
**   o  consider using mode flatten in outputs and modifying groovy functions accordingly
*/

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import org.codehaus.groovy.runtime.StackTraceUtils

/*
** Where to find scripts.
*/
def pipeline_path="/net/gs/vol1/home/bge/eclipse-workspace/bbi-sciatac-analyze"
def script_dir="${pipeline_path}/src"

/*
** Set errorStrategy directive policy.
** Allowed values:
**   "terminate"  (default NextFlow policy)
**   "finish"
**   "ignore"
**   "retry"      (if set to retry, NextFlow retries 1 time (by default: can be changed with 'maxRetries' directive)
**
** Notes:
**   o  the NextFlow documentation describes dynamic retries: errorStrategy { sleep(Math.pow(2, task.attempt) * 200); return 'retry' }
**      this increases exponentially the time between retries, in case of long latencies
*/
def onError = { return( "retry" ) }

/*
** ================================================================================
** Initial set up and checks.
** ================================================================================
*/

/*
** Initial pre-defined, required parameter values.
*/
params.help = false
params.max_cores = 16

/*
** Initialize optional parameters to null.
*/
params.samples = null
params.bowtie_seed = null
params.reads_threshold = null
params.no_secondary = null
params.calculate_banding_scores = null
params.topic_models = null
params.topics = null


/*
** Various internal parameters.
*/

/*
** MakeWindowedGenomeIntervals parameter.
*/
def genomicIntervalWindowSize = 5000

/*
** MakePromoterSumIntervals parameters.
*/
def proximalUpstream = 1000
def proximalDownstream = 500
def peakToTssDistanceThreshold = 30000

/*
** RunAlign parameters.
*/
def max_cores_align = params.max_cores
def total_memory_align = 36 + 0.25 * max_cores_align
def memory_align = total_memory_align / max_cores_align

/*
** makeMergedPeakRegionCountsProcess parameters.
** Notes:
**   o  double check this value
*/
def flanking_distance_merged_peaks_regions = 0

/*
** makeTssRegionCountsProcess parameters.
** Notes:
**   o  double check this value
*/
def flanking_distance_tss_regions = 1000

/*
**getPerBaseCoverageTssProcess parameters.
** Notes:
**   o  double check this value
*/
def flanking_distance_per_base_tss_region = 1000

/*
** Print usage when run with --help parameter.
*/
if (params.help) {
	writeHelp()
	exit( 0 )
}

/*
** Check for required parameters.
*/
assert ( params.demux_dir && params.analyze_dir && params.genomes_json) : "missing config file: use -c CONFIG_FILE.config that includes demux_dir, analyze_dir, and genomes_json"

/*
** Report run parameter values.
*/
reportRunParams( params )

/*
** Check that required directories exist or can be made.
*/
checkDirectories( params )


/*
** ================================================================================
** Read required run information and set up data structures.
** ================================================================================
*/

/*
** Read args.json file in demux_dir.
*/
def argsJson = readArgsJson( params.demux_dir + "/args.json" )

/*
** Get map of sample names (samples in JSON file) where the keys are
** the sample names and the values are the run+lanes with the samples.
** [ '<sample_name>': <list of lanes stored as 'R00n_L00n' strings> ]
** Notes:
**   o  all distinct samples in the sample JSON file
*/
def sampleLaneJsonMap = getSamplesJson( argsJson )

/*
** Check that args.samples, if given, are in json file and
** return a map of sample lanes to process.
*/
def sampleLaneMap = checkArgsSamples( params, sampleLaneJsonMap )

/*
** Make a list of sorted sample names.
*/
def sampleSortedNames = getSortedSampleNames( sampleLaneMap )

/*
** Read genome file paths from json file.
*/
def genomesJson = readGenomesJson( params, argsJson )

/*
** Make a list of names of the genomes required by the samples.
*/
def genomesRequired = findRequiredGenomes( sampleLaneMap, argsJson )

/*
** Make a map of genomes required by the samples.
** [ '<sample_name>': <genome_map_from_json_file> ]
*/
def sampleGenomeMap = getSampleGenomeMap( sampleLaneMap, argsJson, genomesJson )

/*
** ================================================================================
** Additional checks.
** ================================================================================
*/

/*
** Check for trimmed fastq files.
*/
checkFastqs( params, sampleLaneJsonMap )


/*
** Write processing args.json file(s).
** Notes:
**   o  perhaps write sample-specific JSON files
*/
def jsonFilename = params.analyze_dir + '/args.json'
writeRunDataJsonFile( params, argsJson, sampleGenomeMap, jsonFilename )


/*
** ================================================================================
** Prepare output directories and required files.
** ================================================================================
*/

/*
** Sort TSS bed files of required genomes.
*/
Channel
	.fromList( sortTssBedChannelSetup( sampleSortedNames, sampleGenomeMap, genomesJson ) )
	.set { sortTssBedInChannel }

process sortTssBedProcess {
    cache 'lenient'
    errorStrategy onError
	publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "precheck" ) }, pattern: "*.bed.gz", mode: 'copy'
	
	input:
	val tssBedMap from sortTssBedInChannel
	
	output:
	file("${tssBedMap['outBed']}") into sortTssBedOutChannel
	
	script:
	"""
	zcat ${tssBedMap['inBed']} | sort -k1,1V -k2,2n -k3,3n | gzip > ${tssBedMap['outBed']}
	"""
}

sortTssBedOutChannel
    .into { sortTssBedOutChannelCopy01;
            sortTssBedOutChannelCopy02 }

/*
** Sort chromosome size files of required genomes.
*/
Channel
	.fromList( sortChromosomeSizeChannelSetup( sampleSortedNames, sampleGenomeMap, genomesJson ) )
	.set { sortChromosomeSizeInChannel }

process sortChromosomeSizeProcess {
    cache 'lenient'
    errorStrategy errorStrategy { sleep(Math.pow(2, task.attempt) * 200); return 'retry' }
	publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "precheck" ) }, pattern: "*.txt", mode: 'copy'
	
	input:
		val chromosomeSizeMap from sortChromosomeSizeInChannel
	
	output:
		file("${chromosomeSizeMap['outTxt']}") into sortChromosomeSizeOutChannel
	
	script:
	"""
	cat ${chromosomeSizeMap['inTxt']} | sort -k1,1V -k2,2n -k3,3n > ${chromosomeSizeMap['outTxt']}
	"""
}

sortChromosomeSizeOutChannel
    .into { sortChromosomeSizeOutChannelCopy01;
            sortChromosomeSizeOutChannelCopy02;
            sortChromosomeSizeOutChannelCopy03;
            sortChromosomeSizeOutChannelCopy04 }


/*
** ================================================================================
** Run alignments and merge resulting BAM files.
** ================================================================================
*/

/*
** Run alignments.
*/
Channel
	.fromList( runAlignChannelSetup( params, argsJson, sampleLaneMap, genomesJson ) )
	.set { runAlignInChannel }

process runAlignProcess {
	cache 'lenient'
    errorStrategy onError
	penv 'serial'
	cpus max_cores_align
	memory "${memory_align}G"
	module 'java/latest:modules:modules-init:modules-gs:bowtie2/2.2.3:samtools/1.9'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "align_reads" ) }, pattern: "*.bam", mode: 'copy'
 
	input:
		val alignMap from runAlignInChannel

	output:
		file("*.bam") into runAlignOutChannel

	script:
	"""
	bowtie2 -3 1 \
		-X 2000 \
		-p ${max_cores_align} \
		-x ${alignMap['genome_index']} \
		-1 ${alignMap['fastq1']} \
		-2 ${alignMap['fastq2']} ${alignMap['bowtieSeed']} \
		| samtools view -L ${alignMap['whitelist']} -f3 -F12 -q10 -bS - > ${alignMap['bamfile']}.tmp.bam
        samtools sort -T ${alignMap['bamfile']}.sorttemp --threads 4 ${alignMap['bamfile']}.tmp.bam -o ${alignMap['bamfile']}
        rm ${alignMap['bamfile']}.tmp.bam
	"""
}


/*
** Merge alignment bam files.
*/
runAlignOutChannel
	.toList()
	.flatMap { mergeBamChannelSetup( it, sampleLaneMap ) }
	.set { mergeBamsInChannel }

process mergeBamsProcess {
    cache 'lenient'
    errorStrategy onError
	penv 'serial'
	cpus 8
	memory "16G"
	module 'java/latest:modules:modules-init:modules-gs:sambamba/0.6.5:zlib/1.2.6:samtools/1.9'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "merge_bams" ) }, pattern: "*.bam", mode: 'copy'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "merge_bams" ) }, pattern: "*.bai", mode: 'copy'

	input:
		set file(inBams), outBam from mergeBamsInChannel

	output:
		file("*") into mergeBamsOutChannelBam
	
	script:
	"""
	sambamba merge --nthreads 8 ${outBam} ${inBams}
	samtools index ${outBam}
	"""
}


/*
** ================================================================================
** Dedup fragments.
** ================================================================================
*/

mergeBamsOutChannelBam
	.toList()
	.flatMap { getUniqueFragmentsChannelSetup( it, sampleSortedNames ) }
	.set { getUniqueFragmentsInChannel }

process getUniqueFragmentsProcess {
    cache 'lenient'
    errorStrategy onError
	cpus 1
	memory "30G"
	module 'java/latest:modules:modules-init:modules-gs:tabix/0.2.6'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "get_unique_fragments" ) }, pattern: "*-transposition_sites.bed.gz*", mode: 'copy'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "get_unique_fragments" ) }, pattern: "*-fragments.txt.gz*", mode: 'copy'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "get_unique_fragments" ) }, pattern: "*-insert_sizes.txt", mode: 'copy'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "get_unique_fragments" ) }, pattern: "*-duplicate_report.txt", mode: 'copy'
    
	input:
		set file(inBam), file(inBai), uniqueFragmentsMap from getUniqueFragmentsInChannel
		
	output:
		file("*-transposition_sites.bed.gz*") into getUniqueFragmentsOutChannelTranspositionSites  // get both -transposition_sites.bed.gz and -transposition_sites.bed.gz.tbi files
		file("*-fragments.txt.gz*") into getUniqueFragmentOutChannelFragments  // get both -fragments.txt.gz and -fragments.txt.gz.tbi files
		file("*-insert_sizes.txt") into getUniqueFragmentsOutChannelInsertSizeDistribution
		file("*-duplicate_report.txt") into getUniqueFragmentsOutChannelDuplicateReport
		
	script:
	"""
	source ${pipeline_path}/load_python_env_reqs.sh
	source ${script_dir}/python_env/bin/activate
	
	fragments_uncompressed=`echo "${uniqueFragmentsMap['fragments_file']}" | sed 's/.gz//'`
	transposition_sites_uncompressed=`echo "${uniqueFragmentsMap['transposition_sites_file']}" | sed 's/.gz//'`
	python ${script_dir}/get_unique_fragments.py \
		${inBam} \
		--fragments \$fragments_uncompressed \
		--transposition_sites_bed \$transposition_sites_uncompressed \
		--duplicate_read_counts ${uniqueFragmentsMap['duplicate_report']} \
		--insert_sizes ${uniqueFragmentsMap['insert_sizes_file']}

	# Index BAM file / bgzip tabix index for fragments file and transposition_sites BED
	bgzip -f \$fragments_uncompressed
	tabix -p bed ${uniqueFragmentsMap['fragments_file']}

	bgzip -f \$transposition_sites_uncompressed
	tabix -p bed ${uniqueFragmentsMap['transposition_sites_file']}
	"""
}

/*
** Make copies of the transposition_sites_file output
** channel.
*/
getUniqueFragmentsOutChannelTranspositionSites
	.into { getUniqueFragmentsOutChannelTranspositionSitesCopy01;
	        getUniqueFragmentsOutChannelTranspositionSitesCopy02;
	        getUniqueFragmentsOutChannelTranspositionSitesCopy03;
	        getUniqueFragmentsOutChannelTranspositionSitesCopy04;
            getUniqueFragmentsOutChannelTranspositionSitesCopy05;
            getUniqueFragmentsOutChannelTranspositionSitesCopy06;
            getUniqueFragmentsOutChannelTranspositionSitesCopy07 }

	
/*
** ================================================================================
** Call peaks.
** ================================================================================
*/

/*
** Merge alignment bam files.
*/
getUniqueFragmentsOutChannelTranspositionSitesCopy01
	.toList()
	.flatMap { callPeaksChannelSetup( it, sampleLaneMap, sampleGenomeMap ) }
	.set { callPeaksInChannel }

process callPeaksProcess {
    cache 'lenient'
    errorStrategy onError
	memory "16G"
	module 'java/latest:modules:modules-init:modules-gs:python/2.7.3:numpy/1.8.1:setuptools/25.1.1:MACS/2.1.0'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "call_peaks" ) }, pattern: "*-peaks.narrowPeak.gz", mode: 'copy'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "call_peaks" ) }, pattern: "*-peaks.xls", mode: 'copy'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "call_peaks" ) }, pattern: "*-summits.bed", mode: 'copy'

	input:
	set file(inBed), file(inTbi), callPeaksMap from callPeaksInChannel

	output:
	file("*-peaks.narrowPeak.gz") into callPeaksOutChannelNarrowPeak
    file("*-peaks.xls") into callPeaksOutChannelPeaksDeadend1
    file("*-summits.bed") into callPeaksOutChannelDeadend2
	
    script:
	"""
    # We used to add --shift -100 and --extsize, but the regions are now pre-shifted and extended
    # as output by other stages (ajh).
	macs2 callpeak -t ${inBed} \
		-f BED \
		-g ${callPeaksMap['macs_genome']} \
		--nomodel \
		--shift -100 \
		--extsize 200 \
		--keep-dup all \
		--call-summits \
		-n ${callPeaksMap['sample_name']} \
		--outdir ${callPeaksMap['out_dir']} 2> /dev/null

	cat ${callPeaksMap['out_dir']}/${callPeaksMap['sample_name']}_peaks.narrowPeak \
		| sort -k1,1V -k2,2n -k3,3n \
		| cut -f1-3 \
		| gzip > ${callPeaksMap['macs_narrowpeak_file']}

	# rename the files using our convention of <sample_name>-* (macs2 names the file with an underscore after the sample name)
	mv ${callPeaksMap['macs_narrowpeak_file']} ${callPeaksMap['output_narrowpeak_file']}
	mv ${callPeaksMap['macs_xls_file']} ${callPeaksMap['output_xls_file']}
	mv ${callPeaksMap['macs_summits_file']} ${callPeaksMap['output_summits_file']}
	
	rm ${callPeaksMap['out_dir']}/${callPeaksMap['sample_name']}_peaks.narrowPeak
	"""
}

callPeaksOutChannelNarrowPeak
    .into { callPeaksOutChannelNarrowPeakCopy01;
            callPeaksOutChannelNarrowPeakCopy02 }
/*
** ================================================================================
** Set up files for downstream processing.
** ================================================================================
*/

/*
** Merge called peaks.
**
*/
callPeaksOutChannelNarrowPeakCopy01
	.toList()
	.flatMap { makePeakFileChannelSetup( it, sampleSortedNames ) }
	.set { mergePeaksInChannel }


process mergePeaksProcess {
	cache 'lenient'
    errorStrategy onError
	memory "5G"
	module 'java/latest:modules:modules-init:modules-gs:bedtools/2.26.0'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "call_peaks" ) }, pattern: "*-merged_peaks.bed", mode: 'copy'

	input:
		set file(inBed), mergePeaksMap from mergePeaksInChannel
			
	output:
		file("*-merged_peaks.bed") into mergePeaksOutChannel
			
	script:
	"""
    zcat ${inBed} \
        | cut -f1-3 \
        | sort -k1,1V -k2,2n -k3,3n \
        | bedtools merge -i - \
        | sort -k1,1V -k2,2n -k3,3n > ${mergePeaksMap['outBed']}
	"""
}
	
mergePeaksOutChannel
    .into { mergePeaksOutChannelCopy01;
            mergePeaksOutChannelCopy02;
            mergePeaksOutChannelCopy03;
            mergePeaksOutChannelCopy04 }

    
/*
** ================================================================================
** Set up to make site counts.
** ================================================================================
*/

/*
** Make windowed genome intervals.
*/
sortChromosomeSizeOutChannelCopy01
    .toList()
    .flatMap { makeWindowedGenomeIntervalsChannelSetup( it, sampleSortedNames, sampleGenomeMap, genomicIntervalWindowSize ) }
    .set { makeWindowedGenomeIntervalsInChannel }

process makeWindowedGenomeIntervalsProcess {
	cache 'lenient'
    errorStrategy onError
	memory "2G"
	module 'java/latest:modules:modules-init:modules-gs:bedtools/2.26.0'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "make_matrices" ) }, pattern: "*genomic_windows.bed", mode: 'copy'

	input:
        set file( inGenomeSizes ), inGenomeSizesMap from makeWindowedGenomeIntervalsInChannel

	output:
        file("*genomic_windows.bed") into makeWindowedGenomeIntervalsOutChannel
        
	script:
	"""
    bedtools makewindows \
        -g ${inGenomeSizes} \
        -w ${inGenomeSizesMap['windowSize']} \
        > ${inGenomeSizesMap['outGenomicWindows']}
	"""
}


/*
** Make promoter sum intervals.
** Notes:
**   o  the gene intervals is defined using one of two sources (files)
*/
mergePeaksOutChannelCopy01
    .toList()
    .flatMap { makePromoterSumIntervalsChannelSetup( it, sampleSortedNames, sampleGenomeMap, proximalDownstream, proximalUpstream, peakToTssDistanceThreshold ) }
    .set { makePromoterSumIntervalsInChannel }

process makePromoterSumIntervalsProcess {
	cache 'lenient'
    errorStrategy onError
	memory "10G"
	module 'java/latest:modules:modules-init:modules-gs:bedtools/2.26.0'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "make_matrices" ) }, pattern: "*-gene_regions.bed.gz", mode: 'copy'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "make_matrices" ) }, pattern: "*-gene_regions_note.txt", mode: 'copy'

	input:
    set file( inPath ), inMap from makePromoterSumIntervalsInChannel
        
	output:
    file( "*-gene_regions.bed.gz" ) into makePromoterSumIntervalsOutChannel
    file( "*-gene_regions_note.txt" ) into makePromoterSumIntervalsOutChannelDeadend

        script:
        if( inMap['hasGeneScoreBed'] == 1 )
        """
        zcat ${inMap['inBedFile']} | sort -k1,1V -k2,2n -k3,3n | gzip > ${inMap['outBed']}
        echo "gene score bed file was used to define gene regions" > ${inMap['sample']}-gene_regions_note.txt
        """

        else
        """
        bedtools closest \
            -d \
            -a ${inPath} \
            -b <(bedtools window \
                -sw \
                -l ${inMap['proximalUpstream']} \
                -r ${inMap['proximalDownstream']} \
                -a ${inMap['tssFile']} \
                -b ${inPath} \
                | cut -f 1-6 \
                | uniq ) \
        | awk '{{ if (\$10 <= ${inMap['peakToTssDistanceThreshold']}) print \$0 }}' \
        | cut -f 1,2,3,7 \
        | sort -k1,1V -k2,2n -k3,3n \
        | uniq | gzip > ${inMap['outBed']}
        echo "TSS definitions and peak locations was used to define gene regions (check this description)" > ${inMap['sample']}-gene_regions_note.txt
        """
}


/*
** ================================================================================
** Make site counts.
** ================================================================================
*/

/*
** Get peak counts in merged peak regions.
*/
getUniqueFragmentsOutChannelTranspositionSitesCopy02
    .toList()
    .flatMap { makeMergedPeakRegionCountsChannelSetupTranspositionSites( it, sampleSortedNames ) }
    .set { makeMergedPeakRegionCountsInChannelTranspositionSites }
    
mergePeaksOutChannelCopy02
    .toList()
    .flatMap { makeMergedPeakRegionCountsChannelSetupMergedPeaks( it, sampleSortedNames ) }
    .set { makeMergedPeakRegionCountsInChannelMergedPeaks }

sortChromosomeSizeOutChannelCopy02
    .toList()
    .flatMap { makeMergedPeakRegionCountsChannelSetupChromosomeSizes( it, sampleSortedNames, sampleGenomeMap ) }
    .set { makeMergedPeakRegionCountsInChannelChromosomeSizes }
    
process makeMergedPeakRegionCountsProcess {
	cache 'lenient'
    errorStrategy onError
	memory '20 GB'
	module 'java/latest:modules:modules-init:modules-gs:zlib/1.2.6:samtools/1.9:bedtools/2.26.0'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "count_report" ) }, pattern: "*-peak_counts.txt", mode: 'copy'

	input:
		set file( inTranspositionSites ), file( inTbiTranspositionSites ), inTranspositionSitesMap from makeMergedPeakRegionCountsInChannelTranspositionSites
		set file( inMergedPeaks ), inMergedPeaksMap from makeMergedPeakRegionCountsInChannelMergedPeaks
		set file( inChromosomeSizes ), inChromosomeSizesMap from makeMergedPeakRegionCountsInChannelChromosomeSizes

	output:
		file( "*-peak_counts.txt" ) into makeMergedPeakRegionCountsOutChannel

	script:
	"""
	source ${pipeline_path}/load_python_env_reqs.sh
    source ${script_dir}/python_env/bin/activate

    TEMP_REGION_FILE="${inMergedPeaksMap['sample']}-temp_regions.gz"
    
    # TODO simplify this... maybe just have a stage that makes this file for TSS rather than complicating the stage itself
    bedtools slop -i ${inMergedPeaks} -g ${inChromosomeSizes} -b ${flanking_distance_merged_peaks_regions} \
    | bedtools merge -i stdin \
    | gzip > \${TEMP_REGION_FILE}

    python ${script_dir}/get_region_counts.py \
        --transposition_sites_intersect <(bedtools intersect -sorted -a ${inTranspositionSites} -b \${TEMP_REGION_FILE}) \
        --output_file ${inMergedPeaksMap['peakCountsOut']}

    rm \${TEMP_REGION_FILE}
	"""
}


/*
** Get peak counts in TSS regions.
*/
getUniqueFragmentsOutChannelTranspositionSitesCopy03
    .toList()
    .flatMap { makeTssRegionCountsChannelSetupTranspositionSites( it, sampleSortedNames ) }
    .set { makeTssRegionCountsInChannelTranspositionSites }
    
sortTssBedOutChannelCopy01
    .toList()
    .flatMap { makeTssRegionCountsChannelSetupTss( it, sampleSortedNames, sampleGenomeMap ) }
    .set { makeTssRegionCountsInChannelTssRegions }

sortChromosomeSizeOutChannelCopy03
    .toList()
    .flatMap { makeTssRegionCountsChannelSetupChromosomeSizes( it, sampleSortedNames, sampleGenomeMap ) }
    .set { makeTssRegionCountsInChannelChromosomeSizes }

process makeTssRegionCountsProcess {
	cache 'lenient'
    errorStrategy onError
	penv 'serial'
	cpus max_cores_align
	memory "${memory_align}G"
	module 'java/latest:modules:modules-init:modules-gs:zlib/1.2.6:samtools/1.9:bedtools/2.26.0'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "count_report" ) }, pattern: "*-tss_counts.txt", mode: 'copy'

    input:
    set file( inTranspositionSites ), file( inTbiTranspositionSites ), inTranspositionSitesMap from makeTssRegionCountsInChannelTranspositionSites
    set file( inTssRegions ), inTssRegionMap from makeTssRegionCountsInChannelTssRegions
    set file( inChromosomeSizes ), inChromosomeSizesMap from makeTssRegionCountsInChannelChromosomeSizes

    output:
    file( "*-tss_counts.txt" ) into makeTssRegionCountsOutChannel

    script:
    """
    source ${pipeline_path}/load_python_env_reqs.sh
    source ${script_dir}/python_env/bin/activate

    TEMP_REGION_FILE="${inTssRegionMap['sample']}-temp_regions.gz"

    # TODO simplify this... maybe just have a stage that makes this file for TSS rather than complicating the stage itself

    bedtools slop -i ${inTssRegions} -g ${inChromosomeSizes} -b ${flanking_distance_tss_regions} \
    | bedtools merge -i stdin \
    | gzip > \${TEMP_REGION_FILE}

    python ${script_dir}/get_region_counts.py \
        --transposition_sites_intersect <(bedtools intersect -sorted -a ${inTranspositionSites} -b \${TEMP_REGION_FILE}) \
        --output_file ${inTssRegionMap['tssCountsOut']}

    rm \${TEMP_REGION_FILE}
    """
}


/*
** ================================================================================
** Make count reports.
** ================================================================================
*/

/*
** Make count reports.
*/
getUniqueFragmentsOutChannelDuplicateReport
    .toList()
    .flatMap { makeCountReportsChannelSetupDuplicateReport( it, sampleSortedNames ) }
    .set { makeCountReportsInChannelDuplicateReport }
    
makeMergedPeakRegionCountsOutChannel
    .toList()
    .flatMap { makeCountReportsChannelSetupMergedPeakRegionCounts( it, sampleSortedNames ) }
    .set { makeCountReportsInChannelMergedPeakRegionCounts }

makeTssRegionCountsOutChannel
    .toList()
    .flatMap { makeCountReportsChannelSetupTssRegionCounts( it, sampleSortedNames ) }
    .set { makeCountReportsInChannelTssRegionCounts }

process makeCountReportsProcess {
	cache 'lenient'
    errorStrategy onError
	memory "10 GB"
	module 'java/latest:modules:modules-init:modules-gs:R/3.6.1'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "count_report" ) }, pattern: "*-count_report.txt", mode: 'copy'

	input:
    set file( inDuplicateReport ), inDuplicateReportMap from makeCountReportsInChannelDuplicateReport
    set file( inMergedPeakRegionCounts ), inMergedPeakRegionCountsMap from makeCountReportsInChannelMergedPeakRegionCounts
    set file( inTssRegionCounts ), inTssRegionCountsMap from makeCountReportsInChannelTssRegionCounts
    
	output:
	file( "*-count_report.txt" ) into makeCountReportsOutChannel

	script:
	"""
    Rscript ${script_dir}/make_count_report.R ${inDuplicateReport} ${inMergedPeakRegionCounts} ${inTssRegionCounts} ${inDuplicateReportMap['outCountReport']}
	"""
}

makeCountReportsOutChannel
    .into { makeCountReportsOutChannelCopy01;
            makeCountReportsOutChannelCopy02 }


/*
** ================================================================================
** Call cells.
** ================================================================================
*/

/*
** Get cell calls.
*/

makeCountReportsOutChannelCopy01
    .toList()
    .flatMap { callCellsChannelSetup( it, sampleSortedNames ) }
    .set { callCellsInChannel }

process callCellsProcess {
	cache 'lenient'
    errorStrategy onError
	memory "5 GB"
	module 'java/latest:modules:modules-init:modules-gs'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "call_cells" ) }, pattern: "*-called_cells.txt", mode: 'copy'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "call_cells" ) }, pattern: "*-called_cells_whitelist.txt", mode: 'copy'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "call_cells" ) }, pattern: "*-called_cells_stats.json", mode: 'copy'

	input:
	set file( inCountReport ), inCountReportMap from callCellsInChannel

	output:
    file( "*-called_cells.txt" ) into callCellsOutChannelCalledCellsCounts
    file( "*-called_cells_whitelist.txt" ) into callCellsOutChannelCalledCellsWhitelist
    file( "*-called_cells_stats.json" ) into callCellsOutChannelCalledCellsStats

	script:
	"""
    source ${pipeline_path}/load_python_env_reqs.sh
    source ${script_dir}/python_env/bin/activate

    python ${script_dir}/call_cells.py ${inCountReportMap['outCalledCellsCounts']} \
                                       ${inCountReportMap['outCellWhitelist']} \
                                       --fit_metadata ${inCountReportMap['outCallCellsStats']} \
                                       --count_report ${inCountReport} ${inCountReportMap['readsThreshold']}

    mv ${inCountReportMap['outCellWhitelist']} cell_whitelist.txt.tmp
    sort cell_whitelist.txt.tmp > ${inCountReportMap['outCellWhitelist']}
    rm cell_whitelist.txt.tmp
	"""
}

callCellsOutChannelCalledCellsWhitelist
    .into { callCellsOutChannelCalledCellsWhitelistCopy01;
            callCellsOutChannelCalledCellsWhitelistCopy02;
            callCellsOutChannelCalledCellsWhitelistCopy03 }


/*
** ================================================================================
** Get per base coverage in TSS regions.
** ================================================================================
*/

/*
** Get per base coverage TSS region coverage
*/
getUniqueFragmentsOutChannelTranspositionSitesCopy04
    .toList()
    .flatMap { getPerBaseCoverageTssChannelSetupTranspositionSites( it, sampleSortedNames ) }
    .set { getPerBaseCoverageTssInChannelTranspositionSites }
    
sortTssBedOutChannelCopy02
    .toList()
    .flatMap { getPerBaseCoverageTssChannelSetupTss( it, sampleSortedNames, sampleGenomeMap ) }
    .set { getPerBaseCoverageTssInChannelTssRegions }

sortChromosomeSizeOutChannelCopy04
    .toList()
    .flatMap { getPerBaseCoverageTssChannelSetupChromosomeSizes( it, sampleSortedNames, sampleGenomeMap ) }
    .set { getPerBaseCoverageTssInChannelChromosomeSizes }

process getPerBaseCoverageTssProcess {
	cache 'lenient'
    errorStrategy onError
	memory '15 GB'
	module 'java/latest:modules:modules-init:modules-gs:bedtools/2.26.0'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "per_base_tss_region_coverage" ) }, pattern: "*-tss_region_coverage.txt.gz", mode: 'copy'

	input:
    set file( inTranspositionSites ), file( inTbiTranspositionSites ), inTranspositionSitesMap from getPerBaseCoverageTssInChannelTranspositionSites
    set file( inTssRegions ), inTssRegionMap from getPerBaseCoverageTssInChannelTssRegions
    set file( inChromosomeSizes ), inChromosomeSizesMap from getPerBaseCoverageTssInChannelChromosomeSizes

	output:
    file( "*-tss_region_coverage.txt.gz" ) into getPerBaseCoverageTssOutChannel

	script:
	"""
    TEMP_OUTPUT_FILE="${inTssRegionMap['sample']}-temp_file.gz"
    
    # First get 2kb regions surrounding TSSs (not strand-specific here)
    # then calculate per-base coverage with bedtools
    # then write any non-zero entries to a file
    bedtools slop -i ${inTssRegions} -g ${inChromosomeSizes} -b ${flanking_distance_per_base_tss_region}  \
    | bedtools coverage -sorted -d -a stdin -b ${inTranspositionSites} \
    | awk '{{ if (\$8 > 0) print \$0 }}' \
    | gzip > \${TEMP_OUTPUT_FILE}

    # Aggregate per-position coverage over all positions across genes, taking strand into account
    Rscript ${script_dir}/aggregate_per_base_tss_region_counts.R \${TEMP_OUTPUT_FILE} ${inTssRegionMap['outCoverage']}

    rm \${TEMP_OUTPUT_FILE}
	"""
}


/*
** ================================================================================
** Make site count matrices.
** ================================================================================
*/

/*
** Make peak matrix.
*/
getUniqueFragmentsOutChannelTranspositionSitesCopy05
    .toList()
    .flatMap { makePeakMatrixChannelSetupTranspositionSites( it, sampleSortedNames ) }
    .set { makePeakMatrixInChannelTranspositionSites }
    
mergePeaksOutChannelCopy03
    .toList()
    .flatMap { makePeakMatrixChannelSetupMergedPeaks( it, sampleSortedNames ) }
    .set { makePeakMatrixInChannelMergedPeaks }

callCellsOutChannelCalledCellsWhitelistCopy01
    .toList()
    .flatMap { makePeakMatrixChannelSetupCellWhitelist( it, sampleSortedNames ) }
    .set { makePeakMatrixInChannelCellWhitelist }

process makePeakMatrixProcess {
	cache 'lenient'
    errorStrategy onError
	memory "25 GB"
	module 'java/latest:modules:modules-init:modules-gs:bedtools/2.26.0:zlib/1.2.6 pigz/latest'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "make_matrices" ) }, pattern: "*-peak_matrix.mtx.gz", mode: 'copy'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "make_matrices" ) }, pattern: "*.txt", mode: 'copy'

	input:
	set file( inTranspositionSites ), file( inTbiTranspositionSites ), inTranspositionSitesMap from makePeakMatrixInChannelTranspositionSites
	set file( inMergedPeaks ), inMergedPeaksMap from makePeakMatrixInChannelMergedPeaks
	set file( inCellWhitelist ), inCellWhitelistMap from makePeakMatrixInChannelCellWhitelist

	output:
	file( "*-peak_matrix.mtx.gz" ) into makePeakMatrixOutChannel
    file( "*.txt" ) into makePeakMatrixOutChannelDeadend

	script:
	"""
    source ${pipeline_path}/load_python_env_reqs.sh
    source ${script_dir}/python_env/bin/activate

    python ${script_dir}/generate_sparse_matrix.py \
    --transposition_sites_intersect <(bedtools intersect -sorted -a ${inMergedPeaks} -b ${inTranspositionSites} -wa -wb) \
    --intervals ${inMergedPeaks} \
    --cell_whitelist ${inCellWhitelist} \
    --matrix_output ${inMergedPeaksMap['outPeakMatrix']}
	"""
}


/*
** Make window matrix.
*/
getUniqueFragmentsOutChannelTranspositionSitesCopy06
    .toList()
    .flatMap { makeWindowMatrixChannelSetupTranspositionSites( it, sampleSortedNames ) }
    .set { makeWindowMatrixInChannelTranspositionSites }

makeWindowedGenomeIntervalsOutChannel
    .toList()
    .flatMap { makeWindowMatrixChannelSetupWindowedGenomeIntervals( it, sampleSortedNames, sampleGenomeMap ) }
    .set { makeWindowMatrixInChannelWindowedGenomeIntervals }

callCellsOutChannelCalledCellsWhitelistCopy02
    .toList()
    .flatMap { makeWindowMatrixChannelSetupCellWhitelist( it, sampleSortedNames ) }
    .set { makeWindowMatrixInChannelCellWhitelist }

process makeWindowMatrixProcess {
	cache 'lenient'
    errorStrategy onError
	memory "25 GB"
	module 'java/latest:modules:modules-init:modules-gs:bedtools/2.26.0:zlib/1.2.6:pigz/latest'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "make_matrices" ) }, pattern: "*-window_matrix.mtx.gz", mode: 'copy'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "make_matrices" ) }, pattern: "*.txt", mode: 'copy'

	input:
    set file( inTranspositionSites ), file( inTbiTranspositionSites ), inTranspositionSitesMap from makeWindowMatrixInChannelTranspositionSites
    set file( inWindowedIntervals ), inWindowedIntervalsMap from makeWindowMatrixInChannelWindowedGenomeIntervals
    set file( inCellWhitelist ), inCellWhitelistMap from makeWindowMatrixInChannelCellWhitelist

	output:
	file( "*" ) into makeWindowMatrixOutChannel

	script:
	"""
    source ${pipeline_path}/load_python_env_reqs.sh
    source ${script_dir}/python_env/bin/activate

    python ${script_dir}/generate_sparse_matrix.py \
    --transposition_sites_intersect <(bedtools intersect -sorted -a ${inWindowedIntervals} -b ${inTranspositionSites} -wa -wb) \
    --intervals ${inWindowedIntervals} \
    --cell_whitelist ${inCellWhitelist} \
    --matrix_output ${inWindowedIntervalsMap['outWindowMatrix']}
	"""
}


/*
** Make promoter matrix.
*/
getUniqueFragmentsOutChannelTranspositionSitesCopy07
    .toList()
    .flatMap { makePromoterMatrixChannelSetupTranspositionSites( it, sampleSortedNames ) }
    .set { makePromoterMatrixInChannelTranspositionSites }

makePromoterSumIntervalsOutChannel
    .toList()
    .flatMap { makePromoterMatrixChannelSetupGeneRegions( it, sampleSortedNames ) }
    .set { makePromoterMatrixInChannelGeneRegions }

callCellsOutChannelCalledCellsWhitelistCopy03
    .toList()
    .flatMap { makePromoterMatrixChannelSetupCellWhitelist( it, sampleSortedNames ) }
    .set { makePromoterMatrixInChannelCellWhitelist }

process makePromoterMatrixProcess {
	cache 'lenient'
    errorStrategy onError
	memory "25 GB"
	module 'java/latest:modules:modules-init:modules-gs:bedtools/2.26.0:zlib/1.2.6:pigz/latest'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "make_matrices" ) }, pattern: "*-promoter_matrix.mtx.gz", mode: 'copy'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "make_matrices" ) }, pattern: "*.txt", mode: 'copy'

	input:
    set file( inTranspositionSites ), file( inTbiTranspositionSites ), inTranspositionSitesMap from makePromoterMatrixInChannelTranspositionSites
    set file( inGeneRegions ), inGeneRegionsMap from makePromoterMatrixInChannelGeneRegions
    set file( inCellWhitelist ), inCellWhitelistMap from makePromoterMatrixInChannelCellWhitelist

	output:
	file( "*-promoter_matrix.mtx.gz" ) into makePromoterMatrixOutChannel
    file( "*.txt" ) into makePromoterMatrixOutChannelDeadend

	script:
	"""
    source ${pipeline_path}/load_python_env_reqs.sh
    source ${script_dir}/python_env/bin/activate

    python ${script_dir}/generate_sparse_matrix.py \
    --transposition_sites_intersect <(bedtools intersect -sorted -a ${inGeneRegions} -b ${inTranspositionSites} -wa -wb) \
    --intervals ${inGeneRegions} \
    --cell_whitelist ${inCellWhitelist} \
    --matrix_output ${inGeneRegionsMap['outPromoterMatrix']}
	"""
}


/*
** ================================================================================
** Make cell call summary.
** ================================================================================
*/

/*
** Summarize cell calls.
** Andrew's pipeline parameters -> NextFlow channels -> filenames
**
**  count_report                            makeCountReportsOutChannelCopy02                *-count_report.txt"
**  call_cells_stats_files                  callCellsOutChannelCalledCellsStats             *-called_cells_stats.json
**  insert_size_distributions               getUniqueFragmentsOutChannelDuplicateReport     *-insert_sizes.txt
**  peak_call_files                         callPeaksOutChannelNarrowPeakCopy02             *-peaks.narrowPeak.gz
**  merged_peaks                            mergePeaksOutChannelCopy04                      *-merged_peaks.bed
**  per_base_tss_region_coverage_files      getPerBaseCoverageTssOutChannel                 *-tss_region_coverage.txt.gz
**  call_cells_summary_plot                 (this function out channel)                     
**  call_cells_summary_stats                (this function out channel)                     
**  window_matrices                         makeWindowMatrixOutChannel                      *-window_matrix.mtx.gz
**  barnyard                                (set from genome file)                             
*/
makeCountReportsOutChannelCopy02
    .toList()
    .flatMap { summarizeCellCallsSetupCountReports( it, sampleSortedNames ) }
    .set { summarizeCellCallsInChannelCountReports }

callCellsOutChannelCalledCellsStats
    .toList()
    .flatMap { summarizeCellCallsSetupCalledCellStats( it, sampleSortedNames ) }
    .set { summarizeCellCallsInChannelCalledCellStats }

getUniqueFragmentsOutChannelInsertSizeDistribution
    .toList()
    .flatMap { summarizeCellCallsSetupInsertSizeDistribution( it, sampleSortedNames ) }
    .set { summarizeCellCallsInChannelInsertSizeDistribution }
    
callPeaksOutChannelNarrowPeakCopy02
    .toList()
    .flatMap { summarizeCellCallsSetupNarrowPeak( it, sampleSortedNames ) }
    .set { summarizeCellCallsInChannelNarrowPeak }

mergePeaksOutChannelCopy04
    .toList()
    .flatMap { summarizeCellCallsSetupMergedPeaks( it, sampleSortedNames ) }
    .set { summarizeCellCallsInChannelMergedPeaks }
    
getPerBaseCoverageTssOutChannel
    .toList()
    .flatMap { summarizeCellCallsSetupPerBaseCoverageTss( it, sampleSortedNames ) }
    .set { summarizeCellCallsInChannelPerBaseCoverageTss }
    
makeWindowMatrixOutChannel
    .toList()
    .flatMap { summarizeCellCallsSetupWindowMatrix( it, sampleSortedNames ) }
    .set { summarizeCellCallsInChannelWindowMatrix }
    
Channel
    .fromList( summarizeCellCallsSetupTestBarnyard( sampleSortedNames, sampleGenomeMap ) )
    .set { summarizeCellCallsInChannelTestBarnyard }

process summarizeCellCallsProcess {
    cache 'lenient'
    errorStrategy onError
    memory "10 GB"
    module 'java/latest:modules:modules-init:modules-gs'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "summarize_cell_calls" ) }, pattern: "*-called_cells_summary.pdf", mode: 'copy'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "summarize_cell_calls" ) }, pattern: "*-called_cells_summary.stats.txt", mode: 'copy'

    input:
    set file( inCountReports ), inCountReportsMap from summarizeCellCallsInChannelCountReports
    set file( inCalledCellStats ), inCalledCellStatsMap from summarizeCellCallsInChannelCalledCellStats
    set file( inInsertSizes ), inInsertSizesMap from summarizeCellCallsInChannelInsertSizeDistribution
    set file( inNarrowPeaks ), inNarrowPeaksMap from summarizeCellCallsInChannelNarrowPeak
    set file( inMergedPeaks ), inMergedPeaksMap from summarizeCellCallsInChannelMergedPeaks
    set file( inPerBaseCoverageTss ), inPerBaseCoverageTss from summarizeCellCallsInChannelPerBaseCoverageTss
    set file( inWindowMatrix ), file( inWindowMatrixRowNames ), file( inWindowMatrixColNames), inWindowMatrixMap from summarizeCellCallsInChannelWindowMatrix
    val inBarnyardMap from summarizeCellCallsInChannelTestBarnyard
    
    output:
    file( "*-called_cells_summary.pdf" ) into summarizeCellCallsOutChannelCallCellsSummaryPlot
    file( "*-called_cells_summary.stats.txt" ) into summarizeCellCallsOutChannelCallCellsSummaryStats

    script:
    """
    BARNYARD_PARAMS=""
    if [ "${inBarnyardMap['isBarnyard']}" == 1 ]
    then
        BARNYARD_PARAMS="--window_matrices ${inWindowMatrix} --barnyard"
    fi
    
    Rscript ${script_dir}/summarize_cell_calls.R \
        --sample_name ${inCountReportsMap['sample']} \
        --read_count_tables ${inCountReports} \
        --stats_files ${inCalledCellStats} \
        --insert_size_tables ${inInsertSizes} \
        --peak_call_files ${inNarrowPeaks} \
        --merged_peaks ${inMergedPeaks} \
        --per_base_tss_region_coverage_files ${inPerBaseCoverageTss} \
        --plot ${inCountReportsMap['outSummaryPlot']} \
        --output_stats ${inCountReportsMap['outSummaryStats']} \${BARNYARD_PARAMS}
    """
}


/*
** Get per cell insert sizes (optional).
**
process getPerCellInsertSizesProcess {
    cache 'lenient'
    errorStrategy onError
    memory "12 GB"
    module 'java/latest:modules:modules-init:modules-gs:zlib/1.2.6:samtools/1.9'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "xxx" ) }, pattern: "xxx", mode: 'copy'

    input:
        xxx

    output:
        xxx

//        # OPTIONAL BANDING SCORES QC
//        if args.calculate_banding_scores:
//            pipeline.add_job(GetPerCellInsertSizes(fragments_file[i], per_cell_insert_sizes[i], cell_whitelist[i]))
//            pipeline.add_job(GetBandingScores(per_cell_insert_sizes[i], banding_scores[i], cell_whitelist[i]))
//
//class GetPerCellInsertSizes:
//    def __init__(self, fragments_file, output_file, cell_whitelist):
//        self.memory = '12G'
//        self.inputs = [fragments_file, cell_whitelist]
//        self.outputs = [output_file]
//
//        self.command = """
//        module purge
//        module load modules modules-init modules-gs
//        module load zlib/1.2.6
//        module load samtools/1.9
//        source {PIPELINE_PATH}/load_python_env_reqs.sh
//        source {SCRIPTS_DIR}/python_env/bin/activate
//
//        python {SCRIPTS_DIR}/get_insert_size_distribution_per_cell.py {fragments_file} {output_file} --barcodes {cell_whitelist}
//        """.format(PIPELINE_PATH=PIPELINE_PATH, SCRIPTS_DIR=SCRIPTS_DIR, fragments_file=fragments_file, cell_whitelist=cell_whitelist, output_file=output_file)

    script:
    """
    xxx
    """
}
*/


/*
** Get banding scores (optional)..
**
process getBandingScoresProcess {
    cache 'lenient'
    errorStrategy onError
    penv 'serial'
    cpus max_cores_align
    memory "${memory_align}G"
    module 'java/latest:modules:modules-init:modules-gs:xxx'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "xxx" ) }, pattern: "xxx", mode: 'copy'

    input:
        xxx

    output:
        xxx

//        # OPTIONAL BANDING SCORES QC
//        if args.calculate_banding_scores:
//            pipeline.add_job(GetPerCellInsertSizes(fragments_file[i], per_cell_insert_sizes[i], cell_whitelist[i]))
//            pipeline.add_job(GetBandingScores(per_cell_insert_sizes[i], banding_scores[i], cell_whitelist[i]))
//
//class GetBandingScores:
//    def __init__(self, insert_size_file, banding_scores_file, cell_whitelist):
//        self.memory = '12G'
//        self.inputs = [insert_size_file, cell_whitelist]
//        self.outputs = [banding_scores_file]
//
//        self.command = """
//        Rscript {SCRIPTS_DIR}/calculate_nucleosome_banding_scores.R {insert_size_file} {banding_scores_file} --barcodes {cell_whitelist}
//        """.format(SCRIPTS_DIR=SCRIPTS_DIR, insert_size_file=insert_size_file, banding_scores_file=banding_scores_file, cell_whitelist=cell_whitelist)


    script:
    """
    xxx
    """
}
*/


/*
** Call motifs.
**
process callMotifsProcess {
	cache 'lenient'
    errorStrategy onError
	penv 'serial'
	cpus max_cores_align
	memory "${memory_align}G"
	module 'java/latest:modules:modules-init:modules-gs:xxx'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "xxx" ) }, pattern: "xxx", mode: 'copy'

	input:
		xxx

	output:
		xxx

//    # MOTIF CALLING IN PEAKS + MOTIF MATRIX THAT IS PER PEAK SET NOT PER SAMPLE
//    motif_calling_outputs = []
//    run_motif_calling = 'motifs' in GENOME_FILES[args.genome] and 'fasta' in GENOME_FILES[args.genome]
//    if run_motif_calling:
//        for gc_bin in range(0, MOTIF_CALLING_GC_BINS):
//            motifs = GENOME_FILES[args.genome]['motifs']
//            fasta = GENOME_FILES[args.genome]['fasta']
//
//            # Note these are really just temp, so not a huge deal that they aren't declared up top
//            # with everything else
//            output_file = peak_motif_files[gc_bin]
//            motif_calling_outputs.append(output_file)
//
//            pipeline.add_job(PeakMotifs(fasta, merged_peaks, motifs, output_file, gc_bin=gc_bin))
//
//    else:
//        print('The specified genome %s, does not have both a "motifs" and "fasta" entry defined, so skipping motif calls.' % args.genome)
//
//class PeakMotifs:
//    def __init__(self, fasta, peaks, motifs, output_file, gc_bin, pwm_threshold=1e-7):
//        self.memory = '10g'
//        self.inputs = [fasta, peaks, motifs]
//        self.outputs = output_file
//
//        self.command = """
//        module purge
//        module load modules modules-init modules-gs
//        source {PIPELINE_PATH}/load_python_env_reqs.sh
//        source {SCRIPTS_DIR}/python_env/bin/activate
//
//        python {SCRIPTS_DIR}/call_peak_motifs.py {fasta} {peaks} {motifs} {output_file} --gc_bin {gc_bin} --pwm_threshold {pwm_threshold}
//        """.format(PIPELINE_PATH=PIPELINE_PATH, SCRIPTS_DIR=SCRIPTS_DIR, fasta=fasta, peaks=peaks, motifs=motifs, output_file=output_file, gc_bin=gc_bin, pwm_threshold=pwm_threshold)
//
	script:
	"""
	xxx
	"""
}
*/



/*
** Make motif matrix
**
process makeMotifMatrixProcess {
	cache 'lenient'
    errorStrategy onError
	penv 'serial'
	cpus max_cores_align
	memory "${memory_align}G"
	module 'java/latest:modules:modules-init:modules-gs:xxx'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "xxx" ) }, pattern: "xxx", mode: 'copy'

	input:
		xxx

	output:
		xxx

// Note: this appears in a conditional block ' if run_motif_calling:' see above
//        pipeline.add_job(MakeMotifMatrix(motif_calling_outputs,
//                            GENOME_FILES[args.genome]['fasta'],
//                            merged_peaks,
//                            GENOME_FILES[args.genome]['motifs'],
//                            peak_tf_matrix))
//
//
//class MakeMotifMatrix:
//    def __init__(self, peak_motif_files, fasta, peaks, motifs, peak_tf_matrix):
//        self.inputs = peak_motif_files + [fasta, peaks, motifs]
//        self.outputs = [peak_tf_matrix]
//        self.memory = '10g'
//
//        peak_motif_files_string = ' '.join(peak_motif_files)
//
//        self.command = f("""
//        module purge
//        module load modules modules-init modules-gs
//        module load zlib/1.2.6 pigz/latest
//        source {PIPELINE_PATH}/load_python_env_reqs.sh
//        source {SCRIPTS_DIR}/python_env/bin/activate
//
//
//        python {SCRIPTS_DIR}/generate_motif_matrix.py \
//        --peak_motif_files {peak_motif_files_string} \
//        --fasta {fasta} \
//        --peaks {peaks} \
//        --motifs {motifs} \
//        --peak_tf_matrix {peak_tf_matrix}
//        """)

	script:
	"""
	xxx
	"""
}
*/


/*
** Make reduced dimension matrix (optional).
**
process makeReducedDimensionMatrixProcess {
	cache 'lenient'
    errorStrategy onError
	penv 'serial'
	cpus max_cores_align
	memory "${memory_align}G"
	module 'java/latest:modules:modules-init:modules-gs:xxx'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "xxx" ) }, pattern: "xxx", mode: 'copy'

	input:
		xxx

	output:
		xxx

//    if not args.no_secondary:
//        for i,sample in enumerate(all_samples):
//            # REDUCE DIMENSION WITH LSI
//            pipeline.add_job(ReduceDimension(peak_matrices[i],
//                            promoter_matrices[i],
//                            svd_coords[i],
//                            umap_coords[i],
//                            tsne_coords[i],
//                            tfidf_matrices[i],
//                            seurat_objects[i]))
//
//            # CISTOPIC MODELS IF REQUESTED (TIME CONSUMING)
//            if args.topic_models:
//                pipeline.add_job(CisTopicModels(peak_matrices[i], topic_models[i], args.topics))
//
//class ReduceDimension:
//    def __init__(self, peak_matrix, promoter_matrix, svd_coords, umap_coords, tsne_coords, tfidf_matrix, seurat_object=None, svd_dimensions=75, include_svd_1=False, sites_per_cell_threshold=100, remove_top_ntile=0.025):
//        self.memory = '25G'
//        self.inputs = [peak_matrix, promoter_matrix]
//        self.outputs = [svd_coords, umap_coords, tsne_coords, tfidf_matrix]
//
//        if seurat_object is not None:
//            self.outputs = self.outputs + [seurat_object]
//
//        self.command = f("""
//        module load zlib/1.2.6 pigz/latest gcc/8.1.0
//
//        Rscript {SCRIPTS_DIR}/reduce_dimensions.R \
//        {peak_matrix} \
//        {promoter_matrix} \
//        --svd_coords {svd_coords} \
//        --umap_coords {umap_coords} \
//        --tsne_coords {tsne_coords} \
//        --tfidf_matrix {tfidf_matrix} \
//        --svd_dimensions {svd_dimensions} \
//        --sites_per_cell_threshold {sites_per_cell_threshold} \
//        --remove_top_ntile {remove_top_ntile} \
//        --fast_tsne_path {SCRIPTS_DIR}/FIt-SNE/bin/fast_tsne""")
//
//        if seurat_object:
//            self.command = self.command + f(' --seurat_object {seurat_object}')
//        if include_svd_1:
//            self.command = self.command + ' --include_svd_1'


	script:
	"""
	xxx
	"""
}
*/


/*
** Make cis topic models (optional).
**
process makeCisTopicModelsProcess {
	cache 'lenient'
    errorStrategy onError
	penv 'serial'
	cpus max_cores_align
	memory "${memory_align}G"
	module 'java/latest:modules:modules-init:modules-gs:xxx'
    publishDir path: "${params.analyze_dir}", saveAs: { qualifyFilename( it, "xxx" ) }, pattern: "xxx", mode: 'copy'

	input:
		xxx

	output:
		xxx


//    if not args.no_secondary:
//        for i,sample in enumerate(all_samples):
//            # REDUCE DIMENSION WITH LSI
//            pipeline.add_job(ReduceDimension(peak_matrices[i],
//                            promoter_matrices[i],
//                            svd_coords[i],
//                            umap_coords[i],
//                            tsne_coords[i],
//                            tfidf_matrices[i],
//                            seurat_objects[i]))
//
//            # CISTOPIC MODELS IF REQUESTED (TIME CONSUMING)
//            if args.topic_models:
//                pipeline.add_job(CisTopicModels(peak_matrices[i], topic_models[i], args.topics))
//
//class CisTopicModels:
//    def __init__(self, matrix, topic_model, topics):
//        self.memory = '40G'
//        self.inputs = [matrix]
//        self.outputs = [topic_model]
//
//        self.command = """
//        Rscript {SCRIPTS_DIR}/generate_cistopic_model.R {matrix} {topic_model} \
//        --topics {topics} \
//        --site_percent_min 0.01 \
//        --min_cell_nonzero 100
//        """.format(SCRIPTS_DIR=SCRIPTS_DIR, matrix=matrix, topic_model=topic_model, topics=' '.join([str(topic) for topic in topics]))

	script:
	"""
	xxx
	"""
}
*/


/*
** Template
**
process template {
	cache 'lenient'
    errorStrategy onError
	penv 'serial'
	cpus max_cores_align
	memory "${memory_align}G"
	module 'java/latest:modules:modules-init:modules-gs:xxx'

	input:
		xxx

	output:
		xxx

	script:
	"""
	xxx
	"""
}
*/



/*
** ================================================================================
** Start of Groovy support functions.
** ================================================================================
*/


/*
** Given a list of paths, return a map where key=file_name and value=path
** Note: this requires that file names be unique in the list.
*/
def getFileMap( inPaths ) {
    def fileMap = [:]
    inPaths.each { aPath ->
        def aFile = aPath.getFileName().toString()
        fileMap[aFile] = aPath
    }
//    assert inPaths.size() == fileMap.size() : 'duplicate file names are incompatible with this function'
    if( inPaths.size() != fileMap.size() ) {
        println 'Error: getFileMap: duplicate file names are incompatible with this function'
        println '  inPaths.size: ' + inPaths.size()
        println '  fileMap.size: ' + fileMap.size()
        println 'inPaths:'
        inPaths.each { aPath ->
            def aFile = aPath.getFileName().toString()
            println "  ap: aFile: ${aFile}"
            println "  ap: aPath: ${aPath}"
        }
       println 'stack trace:'
        def marker = new Throwable()
        StackTraceUtils.sanitize(marker).stackTrace.eachWithIndex { e, i ->
            println "> $i ${e.toString().padRight(30)} ${e.methodName}"
        }
        System.exit( -1 )
    }

    return( fileMap )
}


/*
** Write help information.
*/
def writeHelp() {
	log.info ''
    log.info 'BBI sci-ATAC-seq Analyzer'
    log.info '--------------------------------'
    log.info ''
    log.info 'For reproducibility, please specify all parameters in a config file'
    log.info 'and running'
    log.info '  nextflow run main.cf -c CONFIG_FILE.config.'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run main.nf -c <CONFIG_FILE_NAME>'
    log.info ''
    log.info 'Help: '
    log.info '    --help                              Show this message and exit.'
    log.info ''
    log.info 'Required parameters (specify in your config file):'
    log.info '    params.demux_dir = DEMUX_DIRECTORY (str)     Path to the demux directory with the demultiplexed fastq files.'
    log.info '    params.analyze_dir = OUTPUT_DIRECTORY (str)    Path to the analysis output directory.'
    log.info '    params.genomes_json = GENOMES_JSON         A json file of genome information for analyses.'
    log.info ''
    log.info 'Optional parameters (specify in your config file):'
    log.info '    params.samples = SAMPLES (str)             A list of samples to include in the analysis. Restrict to a subset of samples in demux output.'
    log.info '    process.bowtie_seed = SEED (int)           Bowtie random number generator seed.'
    log.info '    params.reads_threshold = THRESHOLD (int)   Threshold for how many unique reads to consider a cell. Default is to dynamically pick with mclust.'
	log.info '    params.no_secondary (flag)                 This flag prevents dimensionality reduction and clustering.'
	log.info '    params.calculate_banding_scores (flag)     Add flag to calculate banding scores (take a while to run which is why opt-in here).'
	log.info '    params.topic_models (flag)                 Add flag to generate topic models using cisTopic.'
	log.info '    params.topics (str)                        If generating topic models, this allows the user to specify a list of topic counts to try.'
    log.info '    params.max_cores = 16                      The maximum number of cores to use - fewer will be used if appropriate.'
    log.info '    process.maxForks = 20                      The maximum number of processes to run at the same time on the cluster.'
    log.info '    process.queue = "trapnell-short.q"         The queue on the cluster where the jobs should be submitted.'
    log.info ''
    log.info 'Issues? Contact bge@uw.edu'
}


/*
** Report run parameters
*/
def reportRunParams( params ) {
	String s = ""
	s += String.format( "Run parameters\n" )
	s += String.format( "--------------\n" )
	s += String.format( "Demultiplexed fastq directory:        %s\n", params.demux_dir )
	s += String.format( "Analysis output directory:            %s\n", params.analyze_dir )
	if( params.samples != null ) {
		s += String.format( "Samples to include in analysis:   %s\n", params.samples )
	}
    if( params.bowtie_seed != null ) {
        s += String.format( "Bowtie random number generator seed:  %s\n", params.bowtie_seed )
    }
	if( params.reads_threshold != null ) {
	   s += String.format( "Threshold for # reads/cell:           %s\n", params.reads_threshold )
	}
	if( params.no_secondary != null ) {
	   s += String.format( "Omit dimensionality reduction:        set\n" )
	}
	if( params.calculate_banding_scores != null ) {
	   s += String.format( "Calculate banding scores:             %s\n", params.calculate_banding_scores )
	}
	if( params.topic_models != null ) {
	   s += String.format( "Generate topic models using cisTopic: %s\n", params.topic_models )
	}
    if( params.topics != null ) {
	   s += String.format( "If generating topic models, this allows the user to specify a list of topic counts to try.:   %s\n", params.topics )
	}
	s += String.format( "Maximum cores:                        %d\n", params.max_cores )
	s += String.format( "\n" )
	print( s )
	
	/*
	** Let the operator review and accept or reject the run at this point.
	** This might be irritating to the operator. If so, comment it out.
    */
	def answer = System.console().readLine 'Do you want to continue ([n]/y)? '
	if( answer != 'y' ) {
		System.exit( 0 )
	}
}


/*
** Check directories for existence and/or accessibility.
*/
def checkDirectories( params ) {
	/*
	** Check that demux_dir exists.
	*/
	def dirName = params.demux_dir
	def dirHandle = new File( dirName )
	assert dirHandle.exists() : "unable to find demultiplexed fastq directory $dirName"
	assert dirHandle.canRead() : "unable to read demultiplexed fastq directory $dirName"

	/*
	** Check that either the analyze_dir exists or we can create it.
	*/
	dirName = params.analyze_dir
	dirHandle = new File( dirName )
	if( !dirHandle.exists() ) {
		assert dirHandle.mkdirs() : "unable to create analyze output directory $dirName"
	}
}


/*
** Read args.json file from demux directory.
*/
def readArgsJson( String filePath ) {
	File fileHandle = new File( filePath )
	def jsonSlurper = new JsonSlurper()
	def argsJson = jsonSlurper.parse( fileHandle )
	return( argsJson )
}


/*
** Make a map keyed by sample names and the values are lists of run+lane strings
** that have the sample. For example,
** [ 'sample1': [ 'RUN001_L001', 'RUN001_L002' ],
**   'sample2': [ 'RUN001_L001', 'RUN001_L002' ]
** ]
*/
def getSamplesJson( argsJson ) {
	def i
	def sampleLaneMap = [:]
	def samples
	def runs = argsJson.keySet()
	runs.each { aRun ->
		samples = argsJson[aRun]['samples']
		samples.each { aSample ->
			sampleLaneMap[aSample] = []
		}
	}
	runs.each { aRun ->
		samples = argsJson[aRun]['samples']
		def numLanes = argsJson[aRun]['sequencing_run_metadata']['lanes'].toInteger()
		samples.each { aSample ->
			for( i = 1; i <= numLanes; ++i ) {
				sampleLaneMap[aSample].add( String.format( "%s_L%03d", aRun, i ) )
			}
		}
	}
	samples = sampleLaneMap.keySet()
	samples.each { aSample ->
		sampleLaneMap[aSample].unique()
	}
	return( sampleLaneMap )
}


/*
** Check for expected fastq files.
** Notes:
**   o  each combination of sample and lane has two fastq files
**      with names with the form
**         <sample_name>-RUN<run_number>_L<lane_number>_R<read_number>.trimmed.fastq.gz
**      for example,
**         sample1-RUN001_L001_R1.trimmed.fastq.gz
**         sample1-RUN001_L001_R2.trimmed.fastq.gz
*/
def checkFastqs( params, sampleLaneMap ) {
    def fileName
    def fileHandle
	def dirName = params.demux_dir
	def samples = sampleLaneMap.keySet()
	samples.each { aSample ->
		def lanes = sampleLaneMap[aSample]
		lanes.each { aLane ->
			fileName = dirName + '/' + aSample + '/fastqs_trim/' + String.format( '%s-%s_R1.trimmed.fastq.gz', aSample, aLane )
			fileHandle = new File( fileName )
			assert fileHandle.exists() : "unable to find fasta file ${fileName}"
			assert fileHandle.canRead() : "unable to find fasta file ${fileName}"
			fileName = dirName + '/' + aSample + '/fastqs_trim/' + String.format( '%s-%s_R2.trimmed.fastq.gz', aSample, aLane )
			fileHandle = new File( fileName )
			assert fileHandle.exists() : "unable to find fasta file ${fileName}"
			assert fileHandle.canRead() : "unable to find fasta file ${fileName}"
		}
	}
}


/*
** Write run data JSON file(s).
*/
def writeRunDataJsonFile( params, argsJson, sampleGenomeMap, jsonFilename ) {
    File file_json = new File( jsonFilename )
    analyzeDict = [:]
    analyzeDict['run_date'] = new Date()
    if( params.reads_threshold != null ) {
      analyzeDict['params.reads_threshold'] = params.reads_threshold
    }
    if( params.samples != null ) {
        analyzeDict['params.samples'] = params.samples
    }
    argsJson['ANALYZE'] = analyzeDict
    file_json.write( JsonOutput.prettyPrint( JsonOutput.toJson( argsJson ) ) )
}


/*
** Check for required genomes.
** Notes:
**   o  perhaps check that the required genome files exist...
*/
def readGenomesJson( params, argsJson ) {
	def pathFil = params.genomes_json
	def fileHandle = new File( pathFil )
	assert fileHandle.exists() : "unable to find genomes json file \'${pathFil}\'"
	assert fileHandle.canRead() : "unable to find genomes json file \'${pathFil}\'"
	def jsonSlurper  = new JsonSlurper()
	def genomesJson = jsonSlurper.parse( fileHandle )
	def namGenomesJson = genomesJson.keySet()
	def runs = argsJson.keySet()
	runs.each { aRun ->
		def argsGenomes = argsJson[aRun]['genomes'].values()
		argsGenomes.each { aGenome ->
			assert namGenomesJson.contains( aGenome ) : "genome \'${aGenome}\' not in our genomes json file"
		}
	}
	return( genomesJson )
}


def findRequiredGenomes( sampleLaneMap, argsJson ) {
	def genomesRequired = []
	def runs = argsJson.keySet()
	runs.each { aRun ->
		def samples = argsJson[aRun]['genomes'].keySet()
		samples.each { aSample ->
			genomesRequired.add( argsJson[aRun]['genomes'][aSample] )
		}
	}
	return( genomesRequired.unique() )
}


def getSampleGenomeMap( sampleLaneMap, argsJson, genomesJson ) {
	def sampleGenomeMap = [:]
	def samples = sampleLaneMap.keySet()
	def runs = argsJson.keySet()
	samples.each { aSample ->
		runs.each { aRun ->
			if( argsJson[aRun]['genomes'].containsKey( aSample ) ) {
				def genomeName = argsJson[aRun]['genomes'][aSample]
				sampleGenomeMap[aSample] = genomesJson[genomeName]
			}			
		}
	}
	
	/*
	** diagnostics
	println "getSampleGenomeMap map"
	sampleGenomeMap.keySet().each { aSample ->
		println "  sample: ${aSample}"
		println "  genome_index: ${sampleGenomeMap[aSample]['bowtie_index']}"
		println "  tss: ${sampleGenomeMap[aSample]['tss']}"
	}
	*/
	
	return( sampleGenomeMap )
}


/*
** A closure for 'publishing' a file to a sample-specific sub-directory
*/
def qualifyFilename( aFile, aSubDir ) {
    def aSample = aFile.split( '-' )[0]
    def qualifiedFileName = aSample + '/' + aSubDir + '/' + aFile
    return( qualifiedFileName )
}


/*
** Set up channel to sort TSS bed files.
*/
def sortTssBedChannelSetup( sampleSortedNames, sampleGenomeMap, genomesJson ) {
	def inBed
	def outBed
	def tssBedMaps = []
	sampleSortedNames.each { aSample ->
	   def aGenomeJsonMap = sampleGenomeMap[aSample]
	   def hasTss = aGenomeJsonMap.containsKey( 'tss' )
        if( hasTss ) {
            inBed = aGenomeJsonMap['tss']
            outBed = aSample + '-' + aGenomeJsonMap['name'] + '.tss_file.sorted.bed.gz'
        } else {
            inBed = ''
            outBed = ''
        }
        tssBedMaps.add( [ 'hasTss': hasTss, 'inBed': inBed, 'outBed': outBed ] )	   
	}
	
	/*
	** diagnostic
	println "sortTssBedChannelSetup entries"
	tssBedMaps.each { aSet ->
		println "  hasTss: ${aSet['hasTss']}  inBed: ${aSet['inBed']}  outBed: ${aSet['outBed']}"
	}
    */
	
	return( tssBedMaps )
}


/*
** Set up channel to sort chromosome size files.
*/
def sortChromosomeSizeChannelSetup( sampleSortedNames, sampleGenomeMap, genomesJson ) {
	def inTxt
	def outTxt
	def chromosomeSizeMaps = []
    sampleSortedNames.each { aSample ->
       def aGenomeJsonMap = sampleGenomeMap[aSample]
        def hasChromosomeSizes = aGenomeJsonMap.containsKey( 'chromosome_sizes' )
        if( hasChromosomeSizes ) {
            inTxt = aGenomeJsonMap['chromosome_sizes']
            outTxt = aSample + '-' + aGenomeJsonMap['name'] + '.chromosome_sizes.sorted.txt'
        } else {
            inTxt = ''
            outTxt = ''
        }
        chromosomeSizeMaps.add( [ 'hasChromosomeSizes': hasChromosomeSizes, 'inTxt': inTxt, 'outTxt': outTxt ] )
    }

	/*
	** diagnostics
		println "sortChromosomeSizeChannelSetup entries"
        chromosomeSizeMaps.each { aSet ->
		println "  hasChromosomeSizes: ${aSet['hasChromosomeSizes']}"
		println "  inTxt: ${aSet['inTxt']}"
		println "  outTxt: ${aSet['outTxt']}"
	}
    */
	
	return( chromosomeSizeMaps )
}


/*
** Check that params.samples are in args.json file and
** return a map of sample lanes to process. The map
** looks like
** [ 'sample1': [ 'RUN001_L001', 'RUN001_L002' ],
**   'sample2': [ 'RUN001_L001', 'RUN001_L002' ]
** ]
*/
def checkArgsSamples( params, sampleLaneJsonMap ) {
    def sampleLaneMap = [:]
	if( params.samples != null ) {
		def samplesParams = params.samples
		def samplesJson = sampleLaneJsonMap.keySet()
		samplesParams.each { aSample ->
			assert samplesJson.contains( aSample ) : "sample \'${aSample}\' in params.samples is not in the args.json file"
			sampleLaneMap[aSample].addAll( sampleLaneJsonMap[aSample] )
		}
		samplesParams.each { aSample ->
			sampleLaneMap[aSample].unique()
		}
	} else {
		sampleLaneMap = sampleLaneJsonMap
	}
	return( sampleLaneMap )
}


def getSortedSampleNames( sampleLaneMap ) {
    def sampleNames = sampleLaneMap.keySet()
    return( sampleNames.sort() )
}


/*
** ================================================================================
** Run alignments channel setup functions.
** ================================================================================
*/

/*
** Set up channel to run alignments.
** Notes:
**   o  I do not check that the files exist because I construct the expected files
**      from the sampleLaneMap, which is built from the JSON file and run parameter.
**      If the file is missing, the aligner will stop.
**
** Returns:
**   o  a list of maps. Each map sets up a read alignment process. The map
**      has the key:value pairs
**        key             value description
**        ---             -----------------
**        fastq1          name of fastq file of first reads (read pairs)
**        fastq2          name of fastq file of second reads (read pairs)
**        genome_index    path of aligner index directory/files for the required subject sequence
**        whilelist       path of whitelist region .bed file
**        bamfile         path of output bam file
*/
def runAlignChannelSetup( params, argsJson, sampleLaneMap, genomesJson ) {
	def demuxDir = params.demux_dir
	def bowtieSeed = ''
	if( params.bowtie_seed != null ) {
	   bowtieSeed = "--seed ${params.bowtie_seed}"
	}
	def alignMaps = []
	def samples = sampleLaneMap.keySet()
	samples.each { aSample ->
		def lanes = sampleLaneMap[aSample]
		lanes.each { aLane ->
			def fastq1       = demuxDir + '/' + aSample + '/fastqs_trim/' + String.format( '%s-%s_R1.trimmed.fastq.gz', aSample, aLane )
			def fastq2       = demuxDir + '/' + aSample + '/fastqs_trim/' + String.format( '%s-%s_R2.trimmed.fastq.gz', aSample, aLane )
			def theRun       = aLane.split( '_' )[0]
			def genome       = argsJson[theRun]['genomes'][aSample]
			def genome_index = genomesJson[genome]['bowtie_index']
			def whitelist    = genomesJson[genome]['whitelist_regions']
			def bamfile      = String.format( '%s-%s.bam', aSample, aLane )
			alignMaps.add( [ 'fastq1':fastq1, 'fastq2':fastq2, 'genome_index':genome_index, 'whitelist':whitelist, 'bowtieSeed': bowtieSeed, 'bamfile':bamfile ] )
		}
	}
	
	/*
	** diagnostic
	println "runAlignChannelSetup list"
	alignMaps.each { aMap ->
	  println "fastq1: ${aMap['fastq1']}  fastq2: ${aMap['fastq2']}  genome_index: ${aMap['genome_index']}  whitelist: ${aMap['whitelist']}  bamfile: ${aMap['bamfile']}"
	}
	*/
	
	return( alignMaps )
}


/*
** ================================================================================
** Merge BAM file channel setup functions.
** ================================================================================
*/

/*
** Set up channel to merge bam files.
** Notes:
**   o  input file name format: <sample_name>-<run+lane_id>.bam
**   o  output file name format: <sample_name>.bam
**   o  check for expected files
**   o  return a list of tuples where each tuple consists of an
**      output filename string and a list of input bam file
**      java/NextFlow paths. There is one entry for each sample.
*/
def mergeBamChannelSetup( inPaths, sampleLaneMap ) {
	/*
	** Check for expected BAM files.
	*/
	def filesExpected = []
	def samples = sampleLaneMap.keySet()
	samples.each { aSample ->
		def lanes = sampleLaneMap[aSample]
		lanes.each { aLane ->
			def fileName = aSample + '-' + aLane + '.bam'
			filesExpected.add( fileName )
		}
	}
	def filesFound = []
	inPaths.each { aPath ->
		filesFound.add( aPath.getFileName().toString() )
	}
	filesExpected.each { aFile ->
		assert aFile in filesFound : "missing expected  file \'${aFile}\' in channel"
	}

	/*
	** Gather input bam files.
	** Store them in a map of lists keyed by sample name.
	*/
	def sampleLaneBamMap = [:]
	samples.each { aSample ->
		sampleLaneBamMap[aSample] = []
	}
	inPaths.each { aPath ->
		def sampleName = aPath.getFileName().toString().split( '-' )[0]
		sampleLaneBamMap[sampleName].add( aPath )
	}
	
	/*
	** Set up output channel tuples.
	*/
	def outTuples = []
	samples.each { aSample ->
		def outBam = aSample + '-merged.bam'
		def tuple = new Tuple( sampleLaneBamMap[aSample], outBam )
		outTuples.add( tuple )
	}
	
	/*
	** diagnostic
	println "mergeBamChannelSetup list"
	outTuples.each { aTuple ->
		println "  outBam: ${aTuple[0]}"
		println "  inBams: ${aTuple[1]}"

	}
	*/
		
	return( outTuples )
}


/*
** ================================================================================
** Get unique fragments (dedup) channel setup functions.
** ================================================================================
*/

/*
** Set up channel to get unique alignment fragments; BAM and BAI files(dedup).
** Notes:
** input list of files where file name format: <sample>-merged.bam and <sample>-merged.bam.bai
** output: list of tuples in which each tuple consists of [0] bam filename,
**         [1] bai filename, and [2] map of output filenames
*/
def getUniqueFragmentsChannelSetup( inPaths, sampleSortedNames ) {
	/*
	** Check for expected BAM files.
	*/
	def filesExpected = []
	sampleSortedNames.each { aSample ->
		def fileName = aSample + '-merged.bam'
		filesExpected.add( fileName )
	}
    def filesFound = []
    inPaths.each { aPathPair ->
        aPathPair.each { aPath ->
            filesFound.add( aPath.getFileName().toString() )
        }
    }
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    /*
    ** Deal with a mix of *-merged.bam and *-merged.bam.bai files.
    ** Notes:
    **   o  NextFlow stores the pairs of files as a list of two
    **      files so the elements in inPaths is a list of a list
    **      of paths.
    **   o  I put both of the files in the output channel in order
    **      to be certain that NextFlow makes the required symbolic
    **      links in the work directory.
    */
    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPathPair ->
        aPathPair.each { aPath ->
            def aFile = aPath.getFileName().toString()
            def aSample = aFile.split( '-' )[0]
            if( aFile =~ /merged[.]bam$/ ) {
                pathMap[aSample]['bam'] = aPath
            } else if( aFile =~ /merged[.]bam[.]bai$/ ) {
                pathMap[aSample]['bai'] = aPath
            } else {
                println "Warning: getUniqueFragmentsChannelSetup: unexpected file \'${fileName}\'"
            }
        }
    }
    
	/*
	** Set up output channel tuples.
	*/
    
    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inBam = pathMap[aSample]['bam']
        def inBai = pathMap[aSample]['bai']
        def fragments_file = aSample + '-fragments.txt.gz'
        def transposition_sites_file = aSample + '-transposition_sites.bed.gz'
        def insert_sizes_file = aSample + '-insert_sizes.txt'
        def duplicate_report = aSample + '-duplicate_report.txt'
        def tuple = new Tuple( inBam, inBai, [ 'sample':aSample, 'fragments_file': fragments_file, 'transposition_sites_file': transposition_sites_file, 'insert_sizes_file': insert_sizes_file, 'duplicate_report': duplicate_report ] )
        outTuples.add( tuple )
    }
        
	/*
	** diagnostics

	println "getUniqueFragmentsChannelSetup"
	outTuples.each { tuple ->
		println "  inBam: ${tuple[0]}"
        println "  inBai: ${tuple[1]}"
		println "  fragments_file: ${tuple[2]['fragments_file']}"
		println "  transposition_sites_file: ${tuple[2]['transposition_sites_file']}"
		println "  insert_sizes_file: ${tuple[2]['insert_sizes_file']}"
		println "  duplicate_report: ${tuple[2]['duplicate_report']}"
	}
    */
	
	return( outTuples )
}


/*
** ================================================================================
** Call peaks channel setup functions.
** ================================================================================
*/

/*
** Set up channel to call peaks using MACS2.
*/
def callPeaksChannelSetup( inPaths, sampleLaneMap, sampleGenomeMap ) {
	/*
	** Check for expected bed.gz files.
	*/
	def filesExpected = []
	def samples = sampleLaneMap.keySet()
	samples.each { aSample ->
		def fileName = aSample + '-transposition_sites.bed.gz'
		filesExpected.add( fileName )
	}
    samples.each { aSample ->
        def fileName = aSample + '-transposition_sites.bed.gz.tbi'
        filesExpected.add( fileName )
    }
	def filesFound = []
	inPaths.each { aPathPair ->
		aPathPair.each { aPath ->
			filesFound.add( aPath.getFileName().toString() )
		}
	}
	filesExpected.each { aFile ->
		assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
	}

	/*
	** Deal with a mix of *-transposition_sites.bed.gz and
	** *-transposition_sites.bed.gz.tbi files.
	** Notes:
	**   o  NextFlow stores the pairs of files as a list of two
	**      files so the elements in inPaths is a list of a list
	**      of paths.
	**   o  I put both the -transposition_sites.bed.gz and
	**      -transposition_sites.bed.gz.tbi files in the
	**      output channel in order to be certain that NextFlow
	**      makes the required symbolic links in the work
	**      directory.
	*/
	def pathMap = [:]
	samples.each { aSample ->
		pathMap[aSample] = [:]
	}
	inPaths.each { aPathPair ->
		aPathPair.each { aPath ->
			def aFile = aPath.getFileName().toString()
			def aSample = aFile.split( '-' )[0]
			if( aFile =~ /transposition_sites[.]bed[.]gz$/ ) {
				pathMap[aSample]['bed'] = aPath
			} else if( aFile =~ /transposition_sites[.]bed[.]gz[.]tbi$/ ) {
				pathMap[aSample]['tbi'] = aPath
			} else {
				println "Warning: callPeaksChannelSetup: unexpected file \'${fileName}\'"
			}
		}
	}
	
	def callPeaksTuples = []
	samples.each { aSample ->
		def macs_genome = sampleGenomeMap[aSample]['macs_genome']
		def out_dir = 'call_peaks'
		def macs_narrowpeak_file = "${out_dir}/${aSample}_peaks.narrowPeak.gz"
		def output_narrowpeak_file = "${aSample}-peaks.narrowPeak.gz"
		def macs_xls_file="${out_dir}/${aSample}_peaks.xls"
		def output_xls_file="${aSample}-peaks.xls"
		def macs_summits_file="${out_dir}/${aSample}_summits.bed"
		def output_summits_file="${aSample}-summits.bed"
		def tuple = new Tuple( 
			pathMap[aSample]['bed'],
			pathMap[aSample]['tbi'],
			[ 'sample_name': aSample,
			  'macs_genome': macs_genome,
			  'out_dir': out_dir,
			  'macs_narrowpeak_file': macs_narrowpeak_file,
			  'output_narrowpeak_file': output_narrowpeak_file,
			  'macs_xls_file': macs_xls_file,
			  'output_xls_file': output_xls_file,
			  'macs_summits_file': macs_summits_file,
			  'output_summits_file': output_summits_file
			] )
		callPeaksTuples.add( tuple )
	}
	
	/*
	** diagnostics
	println "callPeaksChannelSetup tuples:"
	callPeaksTuples.each { aTuple ->
		println "  inBed: ${aTuple[0]}"
		println "  inTbi: ${aTuple[1]}"
		println "  sample name: ${aTuple[2]['sample_name']}"
	}
	*/
			
	return( callPeaksTuples )
}


/*
** ================================================================================
** Make peak file channel setup functions.
** ================================================================================
*/

/*
** Set up channel of peak files for downstream analysis.
** Return a list a tuples, a tuple for each sample.
*/
def makePeakFileChannelSetup( inPaths, sampleSortedNames ) {
	/*
	** Check for expected files.
	*/
	def filesExpected = []
	sampleSortedNames.each { aSample ->
        def fileName = aSample + '-peaks.narrowPeak.gz'
        filesExpected.add( fileName )
	}
    def fileMap = getFileMap( inPaths )
	def filesFound = fileMap.keySet()
	filesExpected.each { aFile ->
		assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
	}
	
	def outTuples = []
	sampleSortedNames.each { aSample ->
	   def inBed = aSample + '-peaks.narrowPeak.gz'
	   def outBed = aSample + '-merged_peaks.bed'
	   def tuple = new Tuple( fileMap[inBed], [ 'sample': aSample, 'outBed': outBed ] )
	   outTuples.add( tuple )
	}

	/*
	** diagnostics
	println "makePeakFileChannelSetup file list:"
	outTuples.each { aTuple ->
		def aFile = aTuple[0].getFileName().toString()
        println "  inFile: ${aFile}"
		println "  inPath: ${aTuple[0]}"
		println "  outFile: ${aTuple[1]['outBed']}"
	}
    */

	return( outTuples )
}


/*
** ================================================================================
** Make windowed genome intervals channel setup functions.
** ================================================================================
*/

/*
** Set up channel for making genome intervals for peak counting.
*/
def makeWindowedGenomeIntervalsChannelSetup( inPaths, sampleSortedNames, sampleGenomeMap, windowSize ) {
    /*
    ** Check for expected bed files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def aGenomeJsonMap = sampleGenomeMap[aSample]
        def fileName = aSample + '-' + aGenomeJsonMap['name'] + '.chromosome_sizes.sorted.txt'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def aGenomeJsonMap = sampleGenomeMap[aSample]
        def outGenomicWindows = aSample + '-genomic_windows.bed'
        def inGenomicIntervals = aSample + '-' + aGenomeJsonMap['name'] + '.chromosome_sizes.sorted.txt'
        def tuple = new Tuple( fileMap[inGenomicIntervals], [ 'sample': aSample, 'outGenomicWindows': outGenomicWindows, 'windowSize': windowSize ] )
        outTuples.add( tuple )
    }

    def outMaps = []
    sampleSortedNames.each { aSample ->
        def genomicWindows = aSample + '-genomic_windows.bed'
        def aMap = [:]
        aMap['sample'] = aSample
        aMap['genomicWindows'] = genomicWindows
        aMap['windowSize'] = windowSize
        outMaps.add( aMap )
    }
    
    /*
    ** diagnostics

    println "makeWindowedGenomeIntervalsChannelSetup"
    outTuples.each { aTuple ->
        println "sample: ${aTuple[1]['sample']}"
        println "genomeSizes: ${aTuple[1]['genomeSizes']}"
        println "genomicWindows: ${aTuple[1]['genomicWindows']}"
        println "windowSize: ${aTuple[1]['windowSize']}"
    }
    */
  
    return( outTuples )
}

 
/*
** ================================================================================
** Make promoter sum intervals channel setup functions.
** ================================================================================
*/

/*
** Set up a channel for making gene-related windows for peak counting.
*/
def makePromoterSumIntervalsChannelSetup( inPaths, sampleSortedNames, sampleGenomeMap, proximalDownstream, proximalUpstream, peakToTssDistanceThreshold ) {
    /*
    ** Check for expected bed files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-merged_peaks.bed'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def hasGeneScoreBed
    def inBedPath
    def inBedFile
    def tssFile
    def outBed
    def outTuples = []
    sampleSortedNames.each { aSample ->
        if( sampleGenomeMap[aSample].containsKey( 'gene_score_bed' ) ) {
            hasGeneScoreBed = 1
            inBedPath = 'NA'
            inBedFile = sampleGenomeMap[aSample]['gene_score_bed']
            tssFile = 'NA'
            outBed = aSample + '-gene_regions.bed.gz'
        } else {
            hasGeneScoreBed = 0
            inBedName = aSample + '-merged_peaks.bed'
            inBedPath = fileMap[inBedName]
            inBedFile = 'NA'
            tssFile = sampleGenomeMap[aSample]['tss']
            outBed = aSample + '-gene_regions.bed.gz'            
        }
        def tuple = new Tuple( inBedPath, [ 'sample': aSample,
                                            'hasGeneScoreBed': hasGeneScoreBed,
                                            'tssFile': tssFile,
                                            'inBedFile': inBedFile,
                                            'proximalDownstream': proximalDownstream,
                                            'proximalUpstream': proximalUpstream,
                                            'peakToTssDistanceThreshold': peakToTssDistanceThreshold,
                                            'outBed': outBed ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics

    println 'makePromoterSumIntervalsChannelSetup'
    outTuples.each { aTuple ->
        println "inPath: ${aTuple[0]}"
        println "sample: ${aTuple[1]['sample']}"
    }
    */
    
    return( outTuples )
}


/*
** ================================================================================
** Make merged peak region counts channel setup functions.
** ================================================================================
*/

/*
** Set up a channel for counting sites in merged peaks regions.
** This channel has the transposition sites.
** Notes:
**   o  channel has both <sample_name>-transposition_sites.bed.gz and
**      <sample_name>-transposition_sites.bed.gz.tbi files
**      so the inPaths is a list of lists.
**   o  return a list of lists sorted by sample name
*/
def makeMergedPeakRegionCountsChannelSetupTranspositionSites( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected bed files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-transposition_sites.bed.gz'
        filesExpected.add( fileName )
    }
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-transposition_sites.bed.gz.tbi'
        filesExpected.add( fileName )
    }
    def filesFound = []
    inPaths.each { aPathPair ->
        aPathPair.each { aPath ->
            filesFound.add( aPath.getFileName().toString() )
        }
    }
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPathPair ->
        aPathPair.each { aPath ->
            def aFile = aPath.getFileName().toString()
            def aSample = aFile.split( '-' )[0]
            if( aFile =~ /transposition_sites[.]bed[.]gz$/ ) {
                pathMap[aSample]['bed'] = aPath
            } else if( aFile =~ /transposition_sites[.]bed[.]gz[.]tbi$/ ) {
                pathMap[aSample]['tbi'] = aPath
            } else {
                println "Warning: makeMergedPeakRegionCountsChannelSetupTranspositionSites: unexpected file \'${fileName}\'"
            }
        }
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def tuple = new Tuple( pathMap[aSample]['bed'], pathMap[aSample]['tbi'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    /*
    ** diagnostics
    println "makeMergedPeakRegionCountsChannelSetupTranspositionSites"
    outTuples.each { aTuple ->
        println "inBed: ${aTuple[0]}"
        println "inTbi: ${aTuple[1]}"
    }
    */
    
    return( outTuples )
}


/*
** Set up a channel for counting sites in merged peaks regions.
** This channel has the merged peaks.
** Notes:
**   o  inPath filenames have the form <sample_name>-merged_peaks.bed
**   o  return a list of files sorted by sample name.
*/
def makeMergedPeakRegionCountsChannelSetupMergedPeaks( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected bed files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-merged_peaks.bed'
         filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def peakCountsOut = aSample + '-peak_counts.txt'
        def inMergedPeaks = aSample + '-merged_peaks.bed'
        def tuple = new Tuple( fileMap[inMergedPeaks], [ 'sample': aSample, 'peakCountsOut': peakCountsOut ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "makeMergedPeakRegionCountsChannelSetupMergedPeaks merged peak files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
        
    }
    */
    
    return( outTuples )
}


/*
** Set up a channel for counting sites in merged peaks regions.
** This channel has the chromosome sizes file.
*/
def makeMergedPeakRegionCountsChannelSetupChromosomeSizes( inPaths, sampleSortedNames, sampleGenomeMap ) {
    /*
    ** Check for expected txt files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.chromosome_sizes.sorted.txt'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inFileName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.chromosome_sizes.sorted.txt'
        def tuple = new Tuple( fileMap[inFileName], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    /*
    ** diagnostics
    println "makeMergedPeakRegionCountsChannelSetupChromosomeSizes chromosome sizes files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
    */

    return( outTuples )
}


/*
** ================================================================================
** Make TSS regions counts channel setup functions.
** ================================================================================
*/

/*
** Set up a channel for counting sites in TSS regions.
** This channel has the transposition sites.
** Notes:
**   o  channel has both <sample_name>-transposition_sites.bed.gz and
**      <sample_name>-transposition_sites.bed.gz.tbi files
**      so the inPaths is a list of lists.
**   o  return a list of lists sorted by sample name
*/
def makeTssRegionCountsChannelSetupTranspositionSites( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected bed.gz files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-transposition_sites.bed.gz'
        filesExpected.add( fileName )
    }
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-transposition_sites.bed.gz.tbi'
        filesExpected.add( fileName )
    }
    def filesFound = []
    inPaths.each { aPathPair ->
        aPathPair.each { aPath ->
            filesFound.add( aPath.getFileName().toString() )
        }
    }
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPathPair ->
        aPathPair.each { aPath ->
            def aFile = aPath.getFileName().toString()
            def aSample = aFile.split( '-' )[0]
            if( aFile =~ /transposition_sites[.]bed[.]gz$/ ) {
                pathMap[aSample]['bed'] = aPath
            } else if( aFile =~ /transposition_sites[.]bed[.]gz[.]tbi$/ ) {
                pathMap[aSample]['tbi'] = aPath
            } else {
                println "Warning: makeTssRegionCountsChannelSetupTranspositionSites: unexpected file \'${fileName}\'"
            }
        }
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def tuple = new Tuple( pathMap[aSample]['bed'], pathMap[aSample]['tbi'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    /*
    ** diagnostics
    println "makeTssRegionCountsChannelSetupTranspositionSites"
    outTuples.each { aTuple ->
        println "inBed: ${aTuple[0]}"
        println "inTbi: ${aTuple[1]}"
    }
    */
    
    return( outTuples )
}


/*
** Set up a channel for counting sites in TSS regions.
** This channel has the Tss regions in the genome.
** Notes:
**   o  inPath filenames have the form <sample_name>-<genome_name>.tss_file.sorted.bed.gz
**   o  return a list of files sorted by sample name.
*/
def makeTssRegionCountsChannelSetupTss( inPaths, sampleSortedNames, sampleGenomeMap ) {
    /*
    ** Check for expected bed.gz files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.tss_file.sorted.bed.gz'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def tssCountsOut = aSample + '-tss_counts.txt'
        def inBedName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.tss_file.sorted.bed.gz'
        def tuple = new Tuple( fileMap[inBedName], [ 'sample': aSample, 'tssCountsOut': tssCountsOut ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "makeTssRegionCountsChannelSetupTss TSS region files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
        
    }
    */
    
    return( outTuples )
}


/*
** Set up a channel for counting sites in TSS regions.
** This channel has the chromosome sizes file.
*/
def makeTssRegionCountsChannelSetupChromosomeSizes( inPaths, sampleSortedNames, sampleGenomeMap ) {
    /*
    ** Check for expected txt files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.chromosome_sizes.sorted.txt'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inFileName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.chromosome_sizes.sorted.txt'
        def tuple = new Tuple( fileMap[inFileName], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    /*
    ** diagnostics
    println "makeTssRegionCountsChannelSetupChromosomeSizes chromosome sizes files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
    }
    */

    return( outTuples )
}


/*
** ================================================================================
** Make count reports channel setup functions.
** ================================================================================
*/

/*
** Set up a channel for making the count report.
** This channel contains the duplicate report.
*/
def makeCountReportsChannelSetupDuplicateReport( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected txt files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-duplicate_report.txt'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }
    
    def outTuples = []
    sampleSortedNames.each { aSample ->
       def inTxt = aSample + '-duplicate_report.txt'
       def outCountReport = aSample + '-count_report.txt'
       def tuple = new Tuple( fileMap[inTxt], [ 'sample': aSample, 'outCountReport': outCountReport ] )
       outTuples.add( tuple )
    }

    /*
    ** diagnostics
    println "makeCountReportsChannelSetupDuplicateReport file list:"
    outTuples.each { aTuple ->
        def aFile = aTuple[0].getFileName().toString()
        println "  inFile: ${aFile}"
        println "  inPath: ${aTuple[0]}"
        println "  outFile: ${aTuple[1]['outBed']}"
    }
    */

    return( outTuples )
}


/*
** Set up a channel for making the count report.
** This channel contains the merged peak counts.
*/
def makeCountReportsChannelSetupMergedPeakRegionCounts( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected txt files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-peak_counts.txt'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }
    
    def outTuples = []
    sampleSortedNames.each { aSample ->
       def inTxt = aSample + '-peak_counts.txt'
       def tuple = new Tuple( fileMap[inTxt], [ 'sample': aSample ] )
       outTuples.add( tuple )
    }

    /*
    ** diagnostics
    println "makeCountReportsChannelSetupMergedPeakRegionCounts file list:"
    outTuples.each { aTuple ->
        def aFile = aTuple[0].getFileName().toString()
        println "  inFile: ${aFile}"
        println "  inPath: ${aTuple[0]}"
        println "  outFile: ${aTuple[1]['outBed']}"
    }
    */

    return( outTuples )
}


/*
** Set up a channel for making the count report.
** This channel contains the TSS region counts.
*/
def makeCountReportsChannelSetupTssRegionCounts( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected txt files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-tss_counts.txt'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }
    
    def outTuples = []
    sampleSortedNames.each { aSample ->
       def inTxt = aSample + '-tss_counts.txt'
       def tuple = new Tuple( fileMap[inTxt], [ 'sample': aSample ] )
       outTuples.add( tuple )
    }

    /*
    ** diagnostics
    println "makeCountReportsChannelSetupTssRegionCounts file list:"
    outTuples.each { aTuple ->
        def aFile = aTuple[0].getFileName().toString()
        println "  inFile: ${aFile}"
        println "  inPath: ${aTuple[0]}"
        println "  outFile: ${aTuple[1]['outBed']}"
    }
    */

    return( outTuples )
}


/*
** ================================================================================
** Call cells channel setup functions.
** ================================================================================
*/

/*
** Set up a channel for calling cells.
*/
def callCellsChannelSetup( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected txt files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-count_report.txt'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }
    
    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inTxt = aSample + '-count_report.txt'
        def outCalledCellsCounts = aSample + '-called_cells.txt'
        def outCellWhitelist = aSample + '-called_cells_whitelist.txt'
        def outCallCellsStats = aSample + '-called_cells_stats.json'
        def readsThreshold = ""
        if( params.reads_threshold != null ) {
            readsThreshold = "--reads_threshold ${params.reads_threshold}"
        }
        def tuple = new Tuple( fileMap[inTxt], [ 'sample': aSample,
                                                'readsThreshold': readsThreshold,
                                                'outCalledCellsCounts': outCalledCellsCounts,
                                                'outCellWhitelist': outCellWhitelist,
                                                'outCallCellsStats': outCallCellsStats ] )
        outTuples.add( tuple )
    }

    /*
    ** diagnostics
    println "callCellsChannelSetup file list:"
    outTuples.each { aTuple ->
        def aFile = aTuple[0].getFileName().toString()
        println "  inFile: ${aFile}"
        println "  inPath: ${aTuple[0]}"
        println "  outoutCalledCellsCounts: ${aTuple[1]['outCalledCellsCounts']}"
    }
    */

    return( outTuples )
}


/*
** ================================================================================
** Per base TSS region coverage channel setup functions.
** ================================================================================
*/

/*
** Set up a channel for counting per base coverage.
** This channel has the transposition sites.
** Notes:
**   o  channel has both <sample_name>-transposition_sites.bed.gz and
**      <sample_name>-transposition_sites.bed.gz.tbi files
**      so the inPaths is a list of lists.
**   o  return a list of lists sorted by sample name
*/
def getPerBaseCoverageTssChannelSetupTranspositionSites( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected bed.gz files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-transposition_sites.bed.gz'
        filesExpected.add( fileName )
    }
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-transposition_sites.bed.gz.tbi'
        filesExpected.add( fileName )
    }
    def filesFound = []
    inPaths.each { aPathPair ->
        aPathPair.each { aPath ->
            filesFound.add( aPath.getFileName().toString() )
        }
    }
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPathPair ->
        aPathPair.each { aPath ->
            def aFile = aPath.getFileName().toString()
            def aSample = aFile.split( '-' )[0]
            if( aFile =~ /transposition_sites[.]bed[.]gz$/ ) {
                pathMap[aSample]['bed'] = aPath
            } else if( aFile =~ /transposition_sites[.]bed[.]gz[.]tbi$/ ) {
                pathMap[aSample]['tbi'] = aPath
            } else {
                println "Warning: getPerBaseCoverageTssChannelSetupTranspositionSites: unexpected file \'${fileName}\'"
            }
        }
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def tuple = new Tuple( pathMap[aSample]['bed'], pathMap[aSample]['tbi'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    /*
    ** diagnostics
    println "getPerBaseCoverageTssChannelSetupTranspositionSites"
    outTuples.each { aTuple ->
        println "inBed: ${aTuple[0]}"
        println "inTbi: ${aTuple[1]}"
    }
    */
    
    return( outTuples )
}


/*
** Set up a channel for counting per base sites.
** This channel has the TSS regions.
** Notes:
**   o  inPath filenames have the form <sample_name>-<genome_name>.tss_file.sorted.bed.gz
**   o  return a list of files sorted by sample name.
*/
def getPerBaseCoverageTssChannelSetupTss( inPaths, sampleSortedNames, sampleGenomeMap ) {
    /*
    ** Check for expected bed.gz files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.tss_file.sorted.bed.gz'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def outCoverage = aSample + '-tss_region_coverage.txt.gz'
        def inBedName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.tss_file.sorted.bed.gz'
        def tuple = new Tuple( fileMap[inBedName], [ 'sample': aSample, 'outCoverage': outCoverage ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "getPerBaseCoverageTssChannelSetupTss TSS region files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
        
    }
    */
    
    return( outTuples )
}


/*
** Set up a channel for counting per base sites.
** This channel has the chromosome sizes file.
*/
def getPerBaseCoverageTssChannelSetupChromosomeSizes( inPaths, sampleSortedNames, sampleGenomeMap ) {
    /*
    ** Check for expected txt files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.chromosome_sizes.sorted.txt'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inFileName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.chromosome_sizes.sorted.txt'
        def tuple = new Tuple( fileMap[inFileName], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    /*
    ** diagnostics
    println "getPerBaseCoverageTssChannelSetupChromosomeSizes chromosome sizes files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
    }
    */

    return( outTuples )
}


/*
** ================================================================================
** Make peak matrix channel setup functions.
** ================================================================================
*/

def makePeakMatrixChannelSetupTranspositionSites( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected bed.gz files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-transposition_sites.bed.gz'
        filesExpected.add( fileName )
    }
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-transposition_sites.bed.gz.tbi'
        filesExpected.add( fileName )
    }
    def filesFound = []
    inPaths.each { aPathPair ->
        aPathPair.each { aPath ->
            filesFound.add( aPath.getFileName().toString() )
        }
    }
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPathPair ->
        aPathPair.each { aPath ->
            def aFile = aPath.getFileName().toString()
            def aSample = aFile.split( '-' )[0]
            if( aFile =~ /transposition_sites[.]bed[.]gz$/ ) {
                pathMap[aSample]['bed'] = aPath
            } else if( aFile =~ /transposition_sites[.]bed[.]gz[.]tbi$/ ) {
                pathMap[aSample]['tbi'] = aPath
            } else {
                println "Warning: makePeakMatrixChannelSetupTranspositionSites: unexpected file \'${fileName}\'"
            }
        }
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def tuple = new Tuple( pathMap[aSample]['bed'], pathMap[aSample]['tbi'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    /*
    ** diagnostics
    println "makePeakMatrixChannelSetupTranspositionSites"
    outTuples.each { aTuple ->
        println "inBed: ${aTuple[0]}"
        println "inTbi: ${aTuple[1]}"
    }
    */
    
    return( outTuples )
}


def makePeakMatrixChannelSetupMergedPeaks( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected bed files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-merged_peaks.bed'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def outPeakMatrix = aSample + '-peak_matrix.mtx.gz'
        def inMergedPeaks = aSample + '-merged_peaks.bed'
        def tuple = new Tuple( fileMap[inMergedPeaks], [ 'sample': aSample, 'outPeakMatrix': outPeakMatrix ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "makePeakMatrixChannelSetupMergedPeaks merged peak files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
    }
    */
        
    return( outTuples )
}


def makePeakMatrixChannelSetupCellWhitelist( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected txt files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-called_cells_whitelist.txt'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inCellWhitelist = aSample + '-called_cells_whitelist.txt'
        def tuple = new Tuple( fileMap[inCellWhitelist], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "makePeakMatrixChannelSetupCellWhitelist cell whitelist files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
        
    }
    */

    return( outTuples )
}


/*
** ================================================================================
** Make window matrix channel setup functions.
** ================================================================================
*/

def makeWindowMatrixChannelSetupTranspositionSites( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected bed.gz files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-transposition_sites.bed.gz'
        filesExpected.add( fileName )
    }
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-transposition_sites.bed.gz.tbi'
        filesExpected.add( fileName )
    }
    def filesFound = []
    inPaths.each { aPathPair ->
        aPathPair.each { aPath ->
            filesFound.add( aPath.getFileName().toString() )
        }
    }
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPathPair ->
        aPathPair.each { aPath ->
            def aFile = aPath.getFileName().toString()
            def aSample = aFile.split( '-' )[0]
            if( aFile =~ /transposition_sites[.]bed[.]gz$/ ) {
                pathMap[aSample]['bed'] = aPath
            } else if( aFile =~ /transposition_sites[.]bed[.]gz[.]tbi$/ ) {
                pathMap[aSample]['tbi'] = aPath
            } else {
                println "Warning: makeWindowMatrixChannelSetupTranspositionSite: unexpected file \'${fileName}\'"
            }
        }
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def tuple = new Tuple( pathMap[aSample]['bed'], pathMap[aSample]['tbi'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    /*
    ** diagnostics
    println "makeWindowMatrixChannelSetupTranspositionSite"
    outTuples.each { aTuple ->
        println "inBed: ${aTuple[0]}"
        println "inTbi: ${aTuple[1]}"
    }
    */
    
    return( outTuples )
}


def makeWindowMatrixChannelSetupWindowedGenomeIntervals( inPaths, sampleSortedNames, sampleGenomeMap ) {
    /*
    ** Check for expected bed files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-genomic_windows.bed'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def outWindowMatrix = aSample + '-window_matrix.mtx.gz'
        def inGenomicIntervals = aSample + '-genomic_windows.bed'
        def tuple = new Tuple( fileMap[inGenomicIntervals], [ 'sample': aSample, 'outWindowMatrix': outWindowMatrix ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "makeWindowMatrixChannelSetupWindowedGenomeIntervals windowed genomic intervals files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
        
    }
    */

    return( outTuples )
}


def makeWindowMatrixChannelSetupCellWhitelist( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected txt files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-called_cells_whitelist.txt'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inCellWhitelist = aSample + '-called_cells_whitelist.txt'
        def tuple = new Tuple( fileMap[inCellWhitelist], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "makeWindowMatrixChannelSetupCellWhitelist cell whitelist files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
        
    }
    */

    return( outTuples )
}


/*
** ================================================================================
** Make promoter matrix channel setup functions.
** ================================================================================
*/

def makePromoterMatrixChannelSetupTranspositionSites( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected bed.gz files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-transposition_sites.bed.gz'
        filesExpected.add( fileName )
    }
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-transposition_sites.bed.gz.tbi'
        filesExpected.add( fileName )
    }
    def filesFound = []
    inPaths.each { aPathPair ->
        aPathPair.each { aPath ->
            filesFound.add( aPath.getFileName().toString() )
        }
    }
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPathPair ->
        aPathPair.each { aPath ->
            def aFile = aPath.getFileName().toString()
            def aSample = aFile.split( '-' )[0]
            if( aFile =~ /transposition_sites[.]bed[.]gz$/ ) {
                pathMap[aSample]['bed'] = aPath
            } else if( aFile =~ /transposition_sites[.]bed[.]gz[.]tbi$/ ) {
                pathMap[aSample]['tbi'] = aPath
            } else {
                println "Warning: makePromoterMatrixChannelSetupTranspositionSites: unexpected file \'${fileName}\'"
            }
        }
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def tuple = new Tuple( pathMap[aSample]['bed'], pathMap[aSample]['tbi'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    /*
    ** diagnostics
    println "makePromoterMatrixChannelSetupTranspositionSites"
    outTuples.each { aTuple ->
        println "inBed: ${aTuple[0]}"
        println "inTbi: ${aTuple[1]}"
    }
    */
    
    return( outTuples )
}


def makePromoterMatrixChannelSetupGeneRegions( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected bed files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-gene_regions.bed.gz'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def outPromoterMatrix = aSample + '-promoter_matrix.mtx.gz'
        def inGenomicIntervals = aSample + '-gene_regions.bed.gz'
        def tuple = new Tuple( fileMap[inGenomicIntervals], [ 'sample': aSample, 'outPromoterMatrix': outPromoterMatrix ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "makePromoterMatrixChannelSetupGeneRegions gene regions files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
    }
    */

    return( outTuples )
}


def makePromoterMatrixChannelSetupCellWhitelist( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected txt files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-called_cells_whitelist.txt'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inCellWhitelist = aSample + '-called_cells_whitelist.txt'
        def tuple = new Tuple( fileMap[inCellWhitelist], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "makePromoterMatrixChannelSetupCellWhitelist cell whitelist files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
        
    }
    */

    return( outTuples )
}


/*
** ================================================================================
** Summarize cell calls channel setup functions.
** ================================================================================
*/

def summarizeCellCallsSetupCountReports( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected input paths.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-count_report.txt'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inCountReport = aSample + '-count_report.txt'
        def outSummaryPlot = aSample + '-called_cells_summary.pdf'
        def outSummaryStats = aSample + '-called_cells_summary.stats.txt'
        def tuple = new Tuple( fileMap[inCountReport], [ 'sample': aSample,
                                                         'outSummaryPlot': outSummaryPlot,
                                                         'outSummaryStats': outSummaryStats ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "summarizeCellCallsSetupCountReports files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
        
    }
    */

    return( outTuples )
}


def summarizeCellCallsSetupCalledCellStats( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected input pathss.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-called_cells_stats.json'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inCalledCellStats = aSample + '-called_cells_stats.json'
        def tuple = new Tuple( fileMap[inCalledCellStats], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "summarizeCellCallsSetupCalledCellStats files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
        
    }
    */

    return( outTuples )
}


def summarizeCellCallsSetupInsertSizeDistribution( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected input pathss.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-insert_sizes.txt'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inInsertSizes = aSample + '-insert_sizes.txt'
        def tuple = new Tuple( fileMap[inInsertSizes], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "summarizeCellCallsSetupInsertSizeDistribution files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
        
    }
    */

    return( outTuples )
}


def summarizeCellCallsSetupNarrowPeak( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected input pathss.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-peaks.narrowPeak.gz'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inCalledPeaks = aSample + '-peaks.narrowPeak.gz'
        def tuple = new Tuple( fileMap[inCalledPeaks], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "summarizeCellCallsSetupNarrowPeak files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
        
    }
    */

    return( outTuples )
}


def summarizeCellCallsSetupMergedPeaks( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected input pathss.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-merged_peaks.bed'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inMergedPeaks = aSample + '-merged_peaks.bed'
        def tuple = new Tuple( fileMap[inMergedPeaks], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "summarizeCellCallsSetupMergedPeaks files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
        
    }
    */

    return( outTuples )
}


def summarizeCellCallsSetupPerBaseCoverageTss( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected input pathss.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-tss_region_coverage.txt.gz'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inPerBaseCoverageTss = aSample + '-tss_region_coverage.txt.gz'
        def tuple = new Tuple( fileMap[inPerBaseCoverageTss], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "summarizeCellCallsSetupPerBaseCoverage files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[1]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
        
    }
    */

    return( outTuples )
}


/*
** Set up summarize cell calls channel with window matrix related files.
** Notes:
**   o  there are at least three files for each sample in inPaths.
*/
def summarizeCellCallsSetupWindowMatrix( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected input pathss.
    */
    def filesFound = []
    inPaths.each { aPathSet ->
        aPathSet.each { aPath ->
            filesFound.add( aPath.getFileName().toString() )
        }
    }

    def filesExpected = []
    def fileName
    sampleSortedNames.each { aSample ->
        fileName = aSample + '-window_matrix.mtx.gz'
        filesExpected.add( fileName )
        fileName = aSample + '-window_matrix.rows.txt'
        filesExpected.add( fileName )
        fileName = aSample + '-window_matrix.columns.txt'
        filesExpected.add( fileName )
    }
    filesExpected.each { aFile ->
        assert aFile in filesFound : "missing expected file \'${aFile}\' in channel"
    }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPathSet ->
        aPathSet.each { aPath ->
            def aFile = aPath.getFileName().toString()
            def aSample = aFile.split( '-' )[0]
            if( aFile =~ /window_matrix[.]mtx[.]gz$/ ) {
                pathMap[aSample]['mtx'] = aPath
            } else if( aFile =~ /window_matrix[.]rows[.]txt$/ ) {
                pathMap[aSample]['row'] = aPath
            } else if( aFile =~ /window_matrix[.]columns[.]txt$/ ) {
                pathMap[aSample]['col'] = aPath
            } else {
                println "Warning: makePromoterMatrixChannelSetupTranspositionSites: unexpected file \'${fileName}\'"
            }
        }
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inWindowMatrix = aSample + '-window_matrix.mtx.gz'
        def tuple = new Tuple( pathMap[aSample]['mtx'], pathMap[aSample]['row'], pathMap[aSample]['col'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
    /*
    ** diagnostics
    println "summarizeCellCallsSetupWindowMatrix files"
    outTuples.each { aTuple ->
        def aPath = aTuple[0]
        def aFile = aPath.getFileName().toString()
        println "aSample: ${aTuple[3]['sample']}"
        println "aFile: ${aFile}"
        println "aPath: ${aPath}"
        
    }
    */

    return( outTuples )
}


def summarizeCellCallsSetupTestBarnyard( sampleSortedNames, sampleGenomeMap ) {
    def outMaps = []
    sampleSortedNames.each { aSample ->
        def aMap = [:]
        def isBarnyard = 0
        if( sampleGenomeMap[aSample].containsKey( 'barnyard' ) && sampleGenomeMap[aSample]['barnyard'] == true ) {
            isBarnyard = 1
        }
        aMap = [ 'sample': aSample, 'isBarnyard': isBarnyard ]
        outMaps.add( aMap )
    }
    
    /*
    ** diagnostics
    println "summarizeCellCallsSetupTestBarnyardx files"
    outMaps.keySet().each { aMap ->
        println "aSample: ${aMap['sample']}"
        println "isBarnyard: ${aMap['isBarnyard']}"
    }
    */

    return( outMaps )
}



