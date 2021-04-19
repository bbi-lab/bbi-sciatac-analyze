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
**   o  my effort to establish some definitions that may be
**      fairly consistent with the Nextflow documentation
**
**      term                description
**      ----                -----------
**      process             Nextflow process block
**      process instance    one process block execution (as submitted to the Sun Grid Engine (SGE))
**      data unit           a simple or complex groovy data instance independent of a channel
**      entry item          a simple or complex data unit composing a channel entry
**      channel entry       an item or 'glob' of items submitted to a process instance
**
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
**      causes NextFlow to dump directory/file hashes. It uses
**      hashes when it resumes so one can use this to help diagnose
**      resume problems.
**   o  NextFlow 'consumes' file type content in channels. This means
**      that a channel can be used only once in a script. In order to
**      use a channel's file content in more than one process block, the
**      channel must be duplicated using the .into{} operator. In
**      contrast to the file content, the 'val' content is not
**      consumed so it can be used in more than one process{} block.
**   o  publishDir note: files in the pattern statement must be in an output channel
**   o  take care to escape dollar signs in shell and script (such as awk script)
**      substitutions. These unescaped dollar signs are interpreted as introducing
**      NextFlow/groovy variables in strings. Unfortunately, NextFlow mis-identifies
**      the problematic variable and its location in the code, which complicates
**      greatly finding the error. As an example, which confused me greatly, is
**
**          bedtools slop -i ${inTssRegions} -g ${inChromosomeSizes} -b ${task.ext.params.flanking_distance}  \
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
**     the script string also, which makes sense.) Perhaps one can use '# // ...'?
**  o  regarding files and channels and the resulting lists following the .toList()
**     operator: it appears that when a single instance of a process block puts a
**     single file into an output channel, the file appears as path variable in
**     the list produced by .toList(). So each element in the list returned by
**     .toList() is a single file. However, if the output statement in a process()
**     block places more than one file in the channel, then the files are stored
**     in a list, and the list returned by .toList() is a list of lists of path
**     variables. I suspect that in this case the .flatten() operator (or flatten
**     mode in the output statement) builds a new list in which each element is
**     a path.
**     At this time I do not flatten the channels so I have nested loops to get to
**     the path variables. For the sake of simplifying the Groovy functions that
**     need to deal with the lists of lists, it would be nice to flatten the lists.
**     However, I am concerned about sprinkling flatten statements throughout the
**     code. I feel fairly strongly that I want to not use the flatten mode in the
**     process output definition block because it affects all downstream uses of
**     the channel. (There may be times one wants lists of lists...). One possibility
**     is to use generally the flatten operator immediately before the .toList()
**     operator. I do know what happens if one uses the .flatten() operator on a
**     'flat' channel but it seems reasonable to suppose that it does nothing.
**  o  there is a test/example script that explores channels called channel_01.nf
**     in bbi-sciatac-demux/nextflow
** Notes on debugging Nextflow
**   o  see the function writeChannelEntries() below for debugging channels
**
** Notes on groovy
**   o  scope (these may not be entirely accurate):
**        o  any variable set without a type is a global variable, e.g., 'x = 3'. This
**           can be a source of problems if 'def' is omitted erroneously when defining
**           a list or dictionary in a function that is called more than once.
**        o  'def', or any type, creates a local variable, e.g. 'int x' or 'def x'.
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
**  o  I have tried to keep the comments accurate but I am certain that there
**     are instances in which I cut and pasted and forgot to edit the strings,
**     or I edited the code and some comments/strings no longer describe
**     accurately the resulting code.
**  o  I regret the often long and complicated names; however,
**     I have not noticed a more satisfying alternative. I considered
**     abbreviations but strings of abbreviations in a name may
**     create more ambiguity and difficulty in comprehension. Part
**     of the challenge is distinguishing processes/variables/functions
**     and instances of these that are similar but not identical, of
**     which there are some in this script.
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
**          a description of the file type to the end of the channel
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
**   o  consider using string variables for filenames (define the variables in the process class: try to define processes that create (output) the files)
**   o  update file_map.docx
**   o  write JSON args.json file (perhaps write a sample-specific JSON file to each sample directory) (done)
**   o  fix reads_threshold in call_cells process (done)
**   o  store files in sub-directories (done)
**   o  build for new Shendure processors (no benefit)
**   o  clean repositories (done for now)
**   o  copy processing information to the sample directories
**      so that investigators can describe the processing in
**      a methods section of a paper
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
**   o  check diagnostic strings
**   o  check comments for accuracy
**   o  check conditional code blocks
**   o  test run aggregation
**   o  check makePromoterSumIntervalsProcess
**   o  consider using mode flatten in outputs and modifying groovy functions accordingly
*/

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import org.codehaus.groovy.runtime.StackTraceUtils
import java.nio.file.Files
import java.nio.file.Path 
import java.nio.file.Paths

/*
** Run date/time.
*/
def timeNow = new Date()

/*
** Where to find scripts.
** Note: script_dir needs to be visible within Groovy functions
**       so there is no 'def', which makes it global.
*/
pipeline_path="$workflow.projectDir"
script_dir="${pipeline_path}/src"

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
** Initial pre-defined, required run-specific, command-line parameter values.
*/
params.help = false
params.motif_calling_gc_bins = 25


/*
** Initialize optional parameters to null.
** Notes:
**   boolean values: true/false
*/
params.samples = null
params.bowtie_cpus = 8
params.bowtie_seed = null
params.reads_threshold = null
params.calculate_banding_scores = null
params.doublet_predict = true
params.filter_blacklist_regions = true

/*
** Print usage when run with --help parameter.
*/
if( params.help ) {
	writeHelp()
	exit( 0 )
}

/*
** Check for required parameters.
*/
if( !params.output_dir ) {
	printErr( "Error: missing params.output_dir in parameter file" )
	System.exit( -1 )
}
if( !params.genomes_json ) {
	printErr( "Error: missing params.genomes_json in parameter file" )
	System.exit( -1 )
}

/*
** Define demux_dir and analyze_dir.
*/
output_dir = params.output_dir.replaceAll("/\\z", "")
demux_dir = output_dir + '/demux_out'
analyze_dir = output_dir + '/analyze_out'

/*
** Check that required directories exist or can be made.
*/
checkDirectories( params )

/*
** Report run parameter values.
*/
reportRunParams( params )

/*
** Archive configuration and samplesheet files in demux_dir.
*/
archiveRunFiles( params, timeNow )


/*
** ================================================================================
** Read required run information and set up data structures.
** ================================================================================
*/

/*
** Read args.json file in demux_dir.
*/
def argsJson = readArgsJson( demux_dir + "/args.json" )

/*
** Get map of sample names (samples in JSON file) where the keys are
** the sample names and the values are the run+lanes with the samples.
** [ '<sample_name>': <list of lanes stored as 'R00n_L00n' strings> ]
** Notes:
**   o  all distinct samples in the sample JSON file
*/
def sampleLaneJsonMap = getSamplesJson( argsJson )

/*
** Get a map of merged peak groups keyed by sample name.
*/
def samplePeakGroupMap = getSamplePeakGroupMap( argsJson )

/*
** Get a map of peak files keyed by sample name.
*/
def samplePeakFileMap = getSamplePeakFileMap( argsJson )

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
** Check that peak files exist.
*/
checkPeakFiles( samplePeakFileMap )


/*
** Check for trimmed fastq files.
*/
checkFastqs( params, sampleLaneJsonMap )


/*
** Write processing args.json file(s).
** Notes:
**   o  perhaps write sample-specific JSON files
*/
def jsonFilename = analyze_dir + '/args.json'
writeRunDataJsonFile( params, argsJson, sampleGenomeMap, jsonFilename, timeNow )


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
	publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "precheck" ) }, pattern: "*.tss_file.sorted.bed.gz", mode: 'copy'
//    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "genome_browser" ) }, pattern: "*.bed.gz", mode: 'copy'
	
	input:
	val tssBedMap from sortTssBedInChannel
	
	output:
	file( "*.tss_file.sorted.bed.gz" ) into sortTssBedOutChannel
	
	script:
	"""
	outBed="${tssBedMap['sample']}-${tssBedMap['genome']}.tss_file.sorted.bed.gz"

	zcat ${tssBedMap['inBed']} | sort -k1,1V -k2,2n -k3,3n | gzip > \${outBed}
	"""
}

sortTssBedOutChannel
    .into { sortTssBedOutChannelCopy01;
            sortTssBedOutChannelCopy02;
            sortTssBedOutChannelCopy03 }

/*
** Sort chromosome size files of required genomes.
*/
Channel
	.fromList( sortChromosomeSizeChannelSetup( sampleSortedNames, sampleGenomeMap, genomesJson ) )
	.set { sortChromosomeSizeInChannel }

process sortChromosomeSizeProcess {
    cache 'lenient'
    errorStrategy errorStrategy { sleep(Math.pow(2, task.attempt) * 200); return 'retry' }
	publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "precheck" ) }, pattern: "*.chromosome_sizes.sorted.txt", mode: 'copy'
	
	input:
	val chromosomeSizeMap from sortChromosomeSizeInChannel
	
	output:
	file( "*.chromosome_sizes.sorted.txt" ) into sortChromosomeSizeOutChannel
	
	script:
	"""
	outTxt="${chromosomeSizeMap['sample']}-${chromosomeSizeMap['genome']}.chromosome_sizes.sorted.txt"

	cat ${chromosomeSizeMap['inTxt']} | sort -k1,1V -k2,2n -k3,3n > \${outTxt}
	"""
}

sortChromosomeSizeOutChannel
    .into { sortChromosomeSizeOutChannelCopy01;
            sortChromosomeSizeOutChannelCopy02;
            sortChromosomeSizeOutChannelCopy03;
            sortChromosomeSizeOutChannelCopy04;
            sortChromosomeSizeOutChannelCopy05;
            sortChromosomeSizeOutChannelCopy06;
            sortChromosomeSizeOutChannelCopy07;
            sortChromosomeSizeOutChannelCopy08;
            sortChromosomeSizeOutChannelCopy09 }


/*
** Sort peak files into sample precheck directories.
*/
Channel
    .fromList( peakFileChannelSetup( sampleSortedNames, samplePeakFileMap ) )
    .set { sortPeakFileInChannel }


process sortPeakFileProcess {
    cache 'lenient'
    errorStrategy errorStrategy { sleep(Math.pow(2, task.attempt) * 200); return 'retry' }
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "precheck" ) }, pattern: "*.bed", mode: 'copy'

    input:
    val peakFileMap from sortPeakFileInChannel

    output:
    file( "*.bed" ) into sortPeakFileOutChannel

    script:
    """
    outBed="${peakFileMap['sample']}-${peakFileMap['nameBed']}"
    zcat -f ${peakFileMap['inBed']} | sort -k1,1V -k2,2n -k3,3n > \${outBed}
    """
}


/*
** ================================================================================
** Run alignments and merge resulting BAM files.
** ================================================================================
*/

/*
** Run alignments.
*/
bowtie_cpus = params.bowtie_cpus < 8 ? params.bowtie_cpus : 8
bowtie_memory = ( 36 + 0.25 * bowtie_cpus ) / bowtie_cpus

Channel
	.fromList( runAlignChannelSetup( params, argsJson, sampleLaneMap, genomesJson ) )
	.set { runAlignInChannel }

process runAlignProcess {
    cpus bowtie_cpus
    memory "${bowtie_memory}G"
	cache 'lenient'
    errorStrategy onError
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "align_reads" ) }, pattern: "*.bam", mode: 'copy'
 
	input:
	val alignMap from runAlignInChannel

	output:
	file( "*.bam" ) into runAlignOutChannel

	script:
	"""
	outBam="${alignMap['sample']}-${alignMap['lane']}.bam"

	bowtie2 -3 1 \
		-X 2000 \
		-p ${task.cpus} \
		-x ${alignMap['genome_index']} \
		-1 ${alignMap['fastq1']} \
		-2 ${alignMap['fastq2']} ${alignMap['seed']} \
		| samtools view -L ${alignMap['whitelist']} -f3 -F12 -q10 -bS - > \${outBam}.tmp.bam
        samtools sort -T \${outBam}.sorttemp --threads 4 \${outBam}.tmp.bam -o \${outBam}
        rm \${outBam}.tmp.bam
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
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "merge_bams" ) }, pattern: "*.bam", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "merge_bams" ) }, pattern: "*.bai", mode: 'copy'

	input:
	set file( inBams ), inMergeBamMap from mergeBamsInChannel

	output:
	file( "*" ) into mergeBamsOutChannelBam
	
	script:
	"""
	outBam="${inMergeBamMap['sample']}-merged.bam"
	
	sambamba merge --nthreads 8 \${outBam} ${inBams}
	samtools index \${outBam}

    mkdir -p ${analyze_dir}/${inMergeBamMap['sample']}/genome_browser
    pushd ${analyze_dir}/${inMergeBamMap['sample']}/genome_browser
    ln -sf ../merge_bams/\${outBam} .
    ln -sf ../merge_bams/\${outBam}.bai .
    popd
	"""
}


/*
** ================================================================================
** Dedup fragments.
** ================================================================================
*/

mergeBamsOutChannelBam
    .flatten()
	.toList()
	.flatMap { getUniqueFragmentsChannelSetupBam( it, sampleSortedNames ) }
	.set { getUniqueFragmentsInChannelBam }

sortChromosomeSizeOutChannelCopy08
    .toList()
    .flatMap { getUniqueFragmentsChannelSetupChromosomeSizes( it, sampleSortedNames, sampleGenomeMap ) }
    .set { getUniqueFragmentsInChannelChromosomeSizes }

/*
** Notes: 
**   o  bedGraphToBigWig fails because the chromosome sizes file does not have all sequences used in the alignment, for example, Mt and contigs.
*/
process getUniqueFragmentsProcess {
    cache 'lenient'
    errorStrategy onError
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "get_unique_fragments" ) }, pattern: "*-transposition_sites.bed.gz*", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "get_unique_fragments" ) }, pattern: "*-fragments.txt.gz*", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "get_unique_fragments" ) }, pattern: "*-insert_sizes.txt", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "get_unique_fragments" ) }, pattern: "*-duplicate_report.txt", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "genome_browser" ) },       pattern: "*-transposition_sites.bed.gz*", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "genome_browser" ) },       pattern: "*-read_alignments.bedgraph", mode: 'copy'
//    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "genome_browser" ) },       pattern: "*-read_alignments.bw", mode: 'copy'
    
	input:
	set file( inBam ), file(inBai), inUniqueFragmentsMap from getUniqueFragmentsInChannelBam
    set file( inGenomeSizes ), inGenomeSizesMap from getUniqueFragmentsInChannelChromosomeSizes

	output:
	file( "*-transposition_sites.bed.gz*" ) into getUniqueFragmentsOutChannelTranspositionSites  // get both -transposition_sites.bed.gz and -transposition_sites.bed.gz.tbi files
	file( "*-fragments.txt.gz*" ) into getUniqueFragmentOutChannelFragments  // get both -fragments.txt.gz and -fragments.txt.gz.tbi files
	file( "*-insert_sizes.txt" ) into getUniqueFragmentsOutChannelInsertSizeDistribution
	file( "*-duplicate_report.txt" ) into getUniqueFragmentsOutChannelDuplicateReport
    file( "*-read_alignments.bedgraph" ) into getUniqueFragmentsOutChannelBedGraph
//    file( "*-read_alignments.bw" ) into getUniqueFragmentsOutChannelBigWig
		
	script:
	"""
	source ${pipeline_path}/load_python_env_reqs.sh
	source ${script_dir}/python_env/bin/activate

	outFragments="${inUniqueFragmentsMap['sample']}-fragments.txt"
	outTranspositionSites="${inUniqueFragmentsMap['sample']}-transposition_sites.bed"
	outInsertSizes="${inUniqueFragmentsMap['sample']}-insert_sizes.txt"
	outDuplicateReport="${inUniqueFragmentsMap['sample']}-duplicate_report.txt"
	
	fragments_uncompressed=`echo "${inUniqueFragmentsMap['fragments_file']}" | sed 's/.gz//'`
	transposition_sites_uncompressed=`echo "${inUniqueFragmentsMap['transposition_sites_file']}" | sed 's/.gz//'`

	python ${script_dir}/get_unique_fragments.py \
		${inBam} \
		--fragments \${outFragments} \
		--transposition_sites_bed \${outTranspositionSites} \
		--duplicate_read_counts \${outDuplicateReport} \
		--insert_sizes \${outInsertSizes}

	# Index BAM file / bgzip tabix index for fragments file and transposition_sites BED
	bgzip -f \${outFragments}
	tabix -p bed \${outFragments}.gz

	bgzip -f \${outTranspositionSites}
	tabix -p bed \${outTranspositionSites}.gz

    # Make bedGraph file for genome browser.
    outBedGraph="${inGenomeSizesMap['sample']}-read_alignments.bedgraph"
    outBigWig="${inGenomeSizesMap['sample']}-read_alignments.bw"
    
    bedtools genomecov -ibam ${inBam} -bga | awk 'BEGIN{OFS="\t"}{\$1="chr" \$1}1' > \${outBedGraph}
    # bedGraphToBigWig \${outBedGraph} $inGenomeSizes \${outBigWig}
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
            getUniqueFragmentsOutChannelTranspositionSitesCopy07;
            getUniqueFragmentsOutChannelTranspositionSitesCopy08 }

getUniqueFragmentOutChannelFragments
	.set { getUniqueFragmentOutChannelFragmentsCopy01 }
	
	
/*
** ================================================================================
** Call peaks.
** ================================================================================
*/

getUniqueFragmentsOutChannelTranspositionSitesCopy01
    .flatten()
	.toList()
	.flatMap { callPeaksChannelSetup( it, sampleLaneMap, sampleGenomeMap ) }
	.set { callPeaksInChannel }

process callPeaksProcess {
    cache 'lenient'
    errorStrategy onError
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "call_peaks" ) }, pattern: "*-peaks.narrowPeak.gz", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "call_peaks" ) }, pattern: "*-peaks.xls", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "call_peaks" ) }, pattern: "*-summits.bed", mode: 'copy'

	input:
	set file( inBed ), file(inTbi), inCallPeaksMap from callPeaksInChannel

	output:
	file( "*-peaks.narrowPeak.gz" ) into callPeaksOutChannelNarrowPeak
    file( "*-peaks.xls" ) into callPeaksOutChannelPeaksDeadend1
    file( "*-summits.bed" ) into callPeaksOutChannelDeadend2
	
    script:
	"""
	# MACS uses an underscore between the sample name and the file type
	# so we need to rename output files.
	outDir="call_peaks"
	outMacsNarrowPeak="\${outDir}/${inCallPeaksMap['sample']}_peaks.narrowPeak"
	outNarrowPeak="${inCallPeaksMap['sample']}-peaks.narrowPeak.gz"
	outMacsXls="\${outDir}/${inCallPeaksMap['sample']}_peaks.xls"
	outXls="${inCallPeaksMap['sample']}-peaks.xls"
	outMacsSummits="\${outDir}/${inCallPeaksMap['sample']}_summits.bed"
	outSummits="${inCallPeaksMap['sample']}-summits.bed"
	
    # We used to add --shift -100 and --extsize, but the regions are now pre-shifted and extended
    # as output by other stages (ajh).
	macs2 callpeak -t ${inBed} \
		-f BED \
		-g ${inCallPeaksMap['macs_genome']} \
		--nomodel \
		--shift -100 \
		--extsize 200 \
		--keep-dup all \
		--call-summits \
		-n ${inCallPeaksMap['sample']} \
		--outdir \${outDir} 2> /dev/null

	cat \${outMacsNarrowPeak} \
		| sort -k1,1V -k2,2n -k3,3n \
		| cut -f1-3 \
		| gzip > \${outNarrowPeak}

	mv \${outMacsXls} \${outXls}
	mv \${outMacsSummits} \${outSummits}
	
	rm \${outMacsNarrowPeak}
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
** Notes:
**   o  the mergePeaksByGroupProcess input channel presents
**      sets of bed file paths that define the merged peak
**      groups as specified in the samplesheet.
**   o  merged peak bed files in ${outGroupBed} are copied to
**      sample-specific bed files and the sample-specific files
**      are submitted to the output channel.
**   o  the sample-specific bed files are merged with external
**      bed peak files if the samplesheet gives an external file
**      for the sample. This step is done in the next process
**      block *processMergePeaksByFile*.
*/
callPeaksOutChannelNarrowPeakCopy01
	.toList()
	.flatMap { makePeakByGroupFileChannelSetup( it, sampleSortedNames, samplePeakGroupMap ) }
	.set { mergePeaksByGroupInChannel }


process mergePeaksByGroupProcess {
	cache 'lenient'
    errorStrategy onError
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "call_peaks" ) }, pattern: "*-group_merged_peaks.bed", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "call_peaks" ) }, pattern: "*-merge_by_group_beds.txt", mode: 'copy'
//    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "genome_browser" ) }, pattern: "*-group_merged_peaks.bed", mode: 'copy'

	input:
	set file( 'inBeds' ), mergePeaksMap from mergePeaksByGroupInChannel

	output:
	file( "*-group_merged_peaks.bed" ) into mergePeaksByGroupOutChannel
    file( "*-merge_by_group_beds.txt" ) into mergePeaksByGroupOutChannelBedList

	script:
	"""
	outGroupBed="${mergePeaksMap['group']}-group_merged_peaks_set.bed"

    zcat inBeds* \
        | cut -f1-3 \
        | sort -k1,1V -k2,2n -k3,3n \
        | bedtools merge -i - \
        | sort -k1,1V -k2,2n -k3,3n > \${outGroupBed}

    #
    # Make copies because Nextflow cannot pass on symbolic links.
    #
    for outSample in ${mergePeaksMap['outSampleList']}
    do
      outBed="\${outSample}-group_merged_peaks.bed"
      cp \${outGroupBed} \${outBed}

      outBedList="\${outSample}-merge_by_group_beds.txt"
      echo "${mergePeaksMap['listBeds']}" > \${outBedList}
    done
	"""
}


/*
** The mergePeaksByGroupOutChannel returns paths bundled in lists
** so flatten the channel.
** Notes:
**   o  'mix' in peak files in sortPeakFileOutChannel. The peak
**       files have the sample names prepended so they are
**       distributed correctly for merging.
*/
mergePeaksByGroupOutChannel
    .flatten()
    .mix( sortPeakFileOutChannel )
    .toList()
    .flatMap { makePeakByFileFileChannelSetup( it, sampleSortedNames, samplePeakGroupMap, samplePeakFileMap ) }
    .set { mergePeaksByFileInChannel }


process mergePeaksByFileProcess {
    cache 'lenient'
    errorStrategy onError
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "call_peaks" ) }, pattern: "*-merged_peaks.bed", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "call_peaks" ) }, pattern: "*-merge_by_file_beds.txt", mode: 'copy'
//    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "genome_browser" ) }, pattern: "*-merged_peaks.bed", mode: 'copy'

    input:
    set file( 'inBeds' ), mergePeaksMap from mergePeaksByFileInChannel

    output:
    file( "*-merged_peaks.bed" ) into mergePeaksByFileOutChannel
    file( "*-merge_by_file_beds.txt" ) into mergePeaksByFileOutChannelBedList

    script:
    """
    outBed="${mergePeaksMap['sample']}-merged_peaks.bed"
    cat inBeds* \
        | cut -f1-3 \
        | sort -k1,1V -k2,2n -k3,3n \
        | bedtools merge -i - \
        | sort -k1,1V -k2,2n -k3,3n > \${outBed}

    outBedList="${mergePeaksMap['sample']}-merge_by_file_beds.txt"
    echo "${mergePeaksMap['listBeds']}" > \${outBedList}
    """
}


mergePeaksByFileOutChannel
    .into { mergePeaksByFileOutChannelCopy01;
            mergePeaksByFileOutChannelCopy02;
            mergePeaksByFileOutChannelCopy03;
            mergePeaksByFileOutChannelCopy04;
            mergePeaksByFileOutChannelCopy05;
            mergePeaksByFileOutChannelCopy06;
            mergePeaksByFileOutChannelCopy07 }


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
    .flatMap { makeWindowedGenomeIntervalsChannelSetup( it, sampleSortedNames, sampleGenomeMap ) }
    .set { makeWindowedGenomeIntervalsInChannel }

process makeWindowedGenomeIntervalsProcess {
	cache 'lenient'
    errorStrategy onError
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "make_matrices" ) }, pattern: "*genomic_windows.bed", mode: 'copy'

	input:
    set file( inGenomeSizes ), inGenomeSizesMap from makeWindowedGenomeIntervalsInChannel

	output:
    file( "*genomic_windows.bed" ) into makeWindowedGenomeIntervalsOutChannel
        
	script:
	"""
	outBed="${inGenomeSizesMap['sample']}-genomic_windows.bed"
	
    bedtools makewindows \
        -g ${inGenomeSizes} \
        -w ${task.ext.window_size} \
        > \${outBed}
	"""
}


/*
** Make promoter sum intervals.
** Notes:
**   o  the gene intervals is defined using one of two sources (files)
*/
mergePeaksByFileOutChannelCopy01
    .toList()
    .flatMap { makePromoterSumIntervalsChannelSetup( it, sampleSortedNames, sampleGenomeMap ) }
    .set { makePromoterSumIntervalsInChannel }

process makePromoterSumIntervalsProcess {
	cache 'lenient'
    errorStrategy onError
	module 'openjdk/latest:modules:modules-init:modules-gs:bedtools/2.26.0'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "make_matrices" ) }, pattern: "*-gene_regions.bed.gz", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "make_matrices" ) }, pattern: "*-gene_regions_note.txt", mode: 'copy'

	input:
    set file( inPath ), inMap from makePromoterSumIntervalsInChannel
        
	output:
    file( "*-gene_regions.bed.gz" ) into makePromoterSumIntervalsOutChannel
    file( "*-gene_regions_note.txt" ) into makePromoterSumIntervalsOutChannelDeadend

    script:
    if( inMap['hasGeneScoreBed'] == 1 )
    	"""
    	outBed="${inMap['sample']}-gene_regions.bed.gz"
        
        zcat ${inMap['inBedFile']} | sort -k1,1V -k2,2n -k3,3n | gzip > \${outBed}
        echo "Gene score bed file was used to define gene regions" > ${inMap['sample']}-gene_regions_note.txt
        """
	else
        """
        outBed="${inMap['sample']}-gene_regions.bed.gz"
        
        bedtools closest \
            -d \
            -a ${inPath} \
            -b <(bedtools window \
                -sw \
                -l ${task.ext.proximal_upstream} \
                -r ${task.ext.proximal_downstream} \
                -a ${inMap['tssFile']} \
                -b ${inPath} \
                | cut -f 1-6 \
                | uniq ) \
        | awk '{{ if (\$10 <= ${task.ext.peak_to_tss_distance_threshold} ) print \$0 }}' \
        | cut -f 1,2,3,7 \
        | sort -k1,1V -k2,2n -k3,3n \
        | uniq | gzip > \${outBed}
        echo "TSS definitions and peak locations were used to define gene regions (check this description)" > ${inMap['sample']}-gene_regions_note.txt
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
    .flatten()
    .toList()
    .flatMap { makeMergedPeakRegionCountsChannelSetupTranspositionSites( it, sampleSortedNames ) }
    .set { makeMergedPeakRegionCountsInChannelTranspositionSites }
    
mergePeaksByFileOutChannelCopy02
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
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "count_report" ) }, pattern: "*-peak_counts.txt", mode: 'copy'

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

	outCounts="${inMergedPeaksMap['sample']}-peak_counts.txt"
	tmpRegions="${inMergedPeaksMap['sample']}-temp_regions.gz"
	
    # TODO simplify this... maybe just have a stage that makes this file for TSS rather than complicating the stage itself
    bedtools slop -i ${inMergedPeaks} -g ${inChromosomeSizes} -b ${task.ext.flanking_distance} \
    | bedtools merge -i stdin \
    | gzip > \${tmpRegions}

    python ${script_dir}/get_region_counts.py \
        --transposition_sites_intersect <(bedtools intersect -sorted -a ${inTranspositionSites} -b \${tmpRegions}) \
        --output_file \${outCounts}

    rm \${tmpRegions}
	"""
}


/*
** Get peak counts in TSS regions.
*/
getUniqueFragmentsOutChannelTranspositionSitesCopy03
    .flatten()
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
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "count_report" ) }, pattern: "*-tss_counts.txt", mode: 'copy'

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

	outCounts="${inTssRegionMap['sample']}-tss_counts.txt"
	tmpRegions="${inTssRegionMap['sample']}-temp_regions.gz"
	
    # TODO simplify this... maybe just have a stage that makes this file for TSS rather than complicating the stage itself

    bedtools slop -i ${inTssRegions} -g ${inChromosomeSizes} -b ${task.ext.flanking_distance} \
    | bedtools merge -i stdin \
    | gzip > \${tmpRegions}

    python ${script_dir}/get_region_counts.py \
        --transposition_sites_intersect <(bedtools intersect -sorted -a ${inTranspositionSites} -b \${tmpRegions}) \
        --output_file \${outCounts}

    rm \${tmpRegions}
    """
}


/*
** ================================================================================
** Make count reports.
** ================================================================================
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
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "count_report" ) }, pattern: "*-count_report.txt", mode: 'copy'

	input:
    set file( inDuplicateReport ), inDuplicateReportMap from makeCountReportsInChannelDuplicateReport
    set file( inMergedPeakRegionCounts ), inMergedPeakRegionCountsMap from makeCountReportsInChannelMergedPeakRegionCounts
    set file( inTssRegionCounts ), inTssRegionCountsMap from makeCountReportsInChannelTssRegionCounts
    
	output:
	file( "*-count_report.txt" ) into makeCountReportsOutChannel

	script:
	"""
	outCountReport="${inDuplicateReportMap['sample']}-count_report.txt"
	
    Rscript ${script_dir}/make_count_report.R ${inDuplicateReport} ${inMergedPeakRegionCounts} ${inTssRegionCounts} \${outCountReport}
	"""
}

makeCountReportsOutChannel
    .into { makeCountReportsOutChannelCopy01;
            makeCountReportsOutChannelCopy02;
            makeCountReportsOutChannelCopy03 }


/*
** ================================================================================
** Call cells.
** ================================================================================
*/

if( params.reads_threshold != null ) {
	readsThresholdParameter = "--reads_threshold ${params.reads_threshold}"
}
else {
readsThresholdParameter = ""
}

makeCountReportsOutChannelCopy01
    .toList()
    .flatMap { callCellsChannelSetup( it, sampleSortedNames ) }
    .set { callCellsInChannel }

process callCellsProcess {
	cache 'lenient'
    errorStrategy onError
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "call_cells" ) }, pattern: "*-called_cells.txt", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "call_cells" ) }, pattern: "*-called_cells_whitelist.txt", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "call_cells" ) }, pattern: "*-called_cells_stats.json", mode: 'copy'

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

	outCalledCellsCounts="${inCountReportMap['sample']}-called_cells.txt"
	outCellWhiteList="${inCountReportMap['sample']}-called_cells_whitelist.txt"
	outCallCellsStats="${inCountReportMap['sample']}-called_cells_stats.json"
	
    python ${script_dir}/call_cells.py \${outCalledCellsCounts} \
                                       \${outCellWhiteList} \
                                       --fit_metadata \${outCallCellsStats} \
                                       --count_report ${inCountReport} ${readsThresholdParameter}

    mv \${outCellWhiteList} cell_whitelist.txt.tmp
    sort cell_whitelist.txt.tmp > \${outCellWhiteList}
    rm cell_whitelist.txt.tmp
	"""
}

callCellsOutChannelCalledCellsWhitelist
    .into { callCellsOutChannelCalledCellsWhitelistCopy01;
            callCellsOutChannelCalledCellsWhitelistCopy02;
            callCellsOutChannelCalledCellsWhitelistCopy03;
            callCellsOutChannelCalledCellsWhitelistCopy04 }


/*
** ================================================================================
** Get per base coverage in TSS regions.
** ================================================================================
*/

getUniqueFragmentsOutChannelTranspositionSitesCopy04
    .flatten()
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
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "per_base_tss_region_coverage" ) }, pattern: "*-tss_region_coverage.txt.gz", mode: 'copy'

	input:
    set file( inTranspositionSites ), file( inTbiTranspositionSites ), inTranspositionSitesMap from getPerBaseCoverageTssInChannelTranspositionSites
    set file( inTssRegions ), inTssRegionMap from getPerBaseCoverageTssInChannelTssRegions
    set file( inChromosomeSizes ), inChromosomeSizesMap from getPerBaseCoverageTssInChannelChromosomeSizes

	output:
    file( "*-tss_region_coverage.txt.gz" ) into getPerBaseCoverageTssOutChannel

	script:
	"""
	outCoverage="${inTssRegionMap['sample']}-tss_region_coverage.txt.gz"
	tmpOut="${inTssRegionMap['sample']}-temp_file.gz"
	
    # First get 2kb regions surrounding TSSs (not strand-specific here)
    # then calculate per-base coverage with bedtools
    # then write any non-zero entries to a file
    bedtools slop -i ${inTssRegions} -g ${inChromosomeSizes} -b ${task.ext.flanking_distance}  \
    | bedtools coverage -sorted -d -a stdin -b ${inTranspositionSites} \
    | awk '{{ if (\$8 > 0) print \$0 }}' \
    | gzip > \${tmpOut}

    # Aggregate per-position coverage over all positions across genes, taking strand into account
    Rscript ${script_dir}/aggregate_per_base_tss_region_counts.R \${tmpOut} \${outCoverage}

    rm \${tmpOut}
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
    .flatten()
    .toList()
    .flatMap { makePeakMatrixChannelSetupTranspositionSites( it, sampleSortedNames ) }
    .set { makePeakMatrixInChannelTranspositionSites }
    
mergePeaksByFileOutChannelCopy03
    .toList()
    .flatMap { makePeakMatrixChannelSetupMergedPeaks( it, sampleSortedNames ) }
    .set { makePeakMatrixInChannelMergedPeaks }

callCellsOutChannelCalledCellsWhitelistCopy01
    .toList()
    .flatMap { makePeakMatrixChannelSetupCellWhitelist( it, sampleSortedNames ) }
    .set { makePeakMatrixInChannelCellWhitelist }

sortChromosomeSizeOutChannelCopy05
    .toList()
    .flatMap { makePeakMatrixChannelSetupChromosomeSizes( it, sampleSortedNames, sampleGenomeMap ) }
    .set { makePeakMatrixChannelInChannelChromosomeSizes }

process makePeakMatrixProcess {
	cache 'lenient'
    errorStrategy onError
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "make_matrices" ) }, pattern: "*-peak_matrix.*", mode: 'copy'

	input:
	set file( inTranspositionSites ), file( inTbiTranspositionSites ), inTranspositionSitesMap from makePeakMatrixInChannelTranspositionSites
	set file( inMergedPeaks ), inMergedPeaksMap from makePeakMatrixInChannelMergedPeaks
	set file( inCellWhitelist ), inCellWhitelistMap from makePeakMatrixInChannelCellWhitelist
    set file( inChromosomeSizes ), inChromosomeSizesMap from makePeakMatrixChannelInChannelChromosomeSizes

	output:
	file( "*-peak_matrix.*" ) into makePeakMatrixOutChannel

	script:
	"""
    source ${pipeline_path}/load_python_env_reqs.sh
    source ${script_dir}/python_env/bin/activate

	outPeakMatrix="${inMergedPeaksMap['sample']}-peak_matrix.mtx.gz"
	
    python ${script_dir}/generate_sparse_matrix.py \
    --transposition_sites_intersect <(bedtools intersect -sorted -g ${inChromosomeSizes} -a ${inMergedPeaks} -b ${inTranspositionSites} -wa -wb) \
    --intervals ${inMergedPeaks} \
    --cell_whitelist ${inCellWhitelist} \
    --matrix_output \${outPeakMatrix}
	"""
}


/*
** Make window matrix.
*/
getUniqueFragmentsOutChannelTranspositionSitesCopy06
    .flatten()
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

sortChromosomeSizeOutChannelCopy06
    .toList()
    .flatMap { makeWindowMatrixChannelSetupChromosomeSizes( it, sampleSortedNames, sampleGenomeMap ) }
    .set { makeWindowMatrixChannelInChannelChromosomeSizes }

process makeWindowMatrixProcess {
	cache 'lenient'
    errorStrategy onError
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "make_matrices" ) }, pattern: "*-window_matrix.*", mode: 'copy'

	input:
    set file( inTranspositionSites ), file( inTbiTranspositionSites ), inTranspositionSitesMap from makeWindowMatrixInChannelTranspositionSites
    set file( inWindowedIntervals ), inWindowedIntervalsMap from makeWindowMatrixInChannelWindowedGenomeIntervals
    set file( inCellWhitelist ), inCellWhitelistMap from makeWindowMatrixInChannelCellWhitelist
    set file( inChromosomeSizes ), inChromosomeSizesMap from makeWindowMatrixChannelInChannelChromosomeSizes

	output:
	file( "*-window_matrix.*" ) into makeWindowMatrixOutChannel

	script:
	"""
    source ${pipeline_path}/load_python_env_reqs.sh
    source ${script_dir}/python_env/bin/activate

	outWindowMatrix="${inWindowedIntervalsMap['sample']}-window_matrix.mtx.gz"
	
    python ${script_dir}/generate_sparse_matrix.py \
    --transposition_sites_intersect <(bedtools intersect -sorted -g ${inChromosomeSizes} -a ${inWindowedIntervals} -b ${inTranspositionSites} -wa -wb) \
    --intervals ${inWindowedIntervals} \
    --cell_whitelist ${inCellWhitelist} \
    --matrix_output \${outWindowMatrix}
	"""
}


/*
** Make promoter matrix.
*/
getUniqueFragmentsOutChannelTranspositionSitesCopy07
    .flatten()
    .toList()
    .flatMap { makePromoterMatrixChannelSetupTranspositionSites( it, sampleSortedNames ) }
    .set { makePromoterMatrixInChannelTranspositionSites }

makePromoterSumIntervalsOutChannel
    .toList()
    .flatMap { makePromoterMatrixChannelSetupGeneRegions( it, sampleSortedNames, sampleGenomeMap ) }
    .set { makePromoterMatrixInChannelGeneRegions }

callCellsOutChannelCalledCellsWhitelistCopy03
    .toList()
    .flatMap { makePromoterMatrixChannelSetupCellWhitelist( it, sampleSortedNames ) }
    .set { makePromoterMatrixInChannelCellWhitelist }

sortChromosomeSizeOutChannelCopy07
    .toList()
    .flatMap { makePromoterMatrixChannelSetupChromosomeSizes( it, sampleSortedNames, sampleGenomeMap ) }
    .set { makePromoterMatrixChannelInChannelChromosomeSizes }

process makePromoterMatrixProcess {
	cache 'lenient'
    errorStrategy onError
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "make_matrices" ) }, pattern: "*-promoter_matrix.*", mode: 'copy'

	input:
    set file( inTranspositionSites ), file( inTbiTranspositionSites ), inTranspositionSitesMap from makePromoterMatrixInChannelTranspositionSites
    set file( inGeneRegions ), inGeneRegionsMap from makePromoterMatrixInChannelGeneRegions
    set file( inCellWhitelist ), inCellWhitelistMap from makePromoterMatrixInChannelCellWhitelist
    set file( inChromosomeSizes ), inChromosomeSizesMap from makePromoterMatrixChannelInChannelChromosomeSizes
    
	output:
	file( "*-promoter_matrix.*" ) into makePromoterMatrixOutChannel

	script:
	"""
    source ${pipeline_path}/load_python_env_reqs.sh
    source ${script_dir}/python_env/bin/activate

	outPromoterMatrix="${inGeneRegionsMap['sample']}-promoter_matrix.mtx.gz"
	
    python ${script_dir}/generate_sparse_matrix.py \
    --transposition_sites_intersect <(bedtools intersect -sorted -g ${inChromosomeSizes} -a ${inGeneRegions} -b ${inTranspositionSites} -wa -wb) \
    --intervals ${inGeneRegions} \
    --cell_whitelist ${inCellWhitelist} \
    --matrix_output \${outPromoterMatrix}
    
    mv ${inGeneRegionsMap['inPromoterMatrixRows']} promoter_matrix_rows.txt.no_metadata
    ${script_dir}/add_gene_metadata.py --in_promoter_matrix_row_name_file promoter_matrix_rows.txt.no_metadata \
                                       --gene_metadata_file ${inGeneRegionsMap['inGeneBodiesGeneMap']} \
                                       --out_promoter_matrix_row_name_file ${inGeneRegionsMap['inPromoterMatrixRows']}
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
**  merged_peaks                            mergePeaksByFileOutChannelCopy04                *-merged_peaks.bed
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

mergePeaksByFileOutChannelCopy04
    .toList()
    .flatMap { summarizeCellCallsSetupMergedPeaks( it, sampleSortedNames ) }
    .set { summarizeCellCallsInChannelMergedPeaks }
    
getPerBaseCoverageTssOutChannel
    .toList()
    .flatMap { summarizeCellCallsSetupPerBaseCoverageTss( it, sampleSortedNames ) }
    .set { summarizeCellCallsInChannelPerBaseCoverageTss }
    
makeWindowMatrixOutChannel
    .flatten()
    .toList()
    .flatMap { summarizeCellCallsSetupWindowMatrix( it, sampleSortedNames ) }
    .set { summarizeCellCallsInChannelWindowMatrix }
    
Channel
    .fromList( summarizeCellCallsSetupTestBarnyard( sampleSortedNames, sampleGenomeMap ) )
    .set { summarizeCellCallsInChannelTestBarnyard }

process summarizeCellCallsProcess {
    cache 'lenient'
    errorStrategy onError
    module 'openjdk/latest:modules:modules-init:modules-gs'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "summarize_cell_calls" ) }, pattern: "*-called_cells_summary.pdf", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "summarize_cell_calls" ) }, pattern: "*-called_cells_summary.stats.txt", mode: 'copy'
    publishDir path: "${output_dir}/analyze_dash/img", pattern: "*.png", mode: 'copy'

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
    file( "*.png") into summarizeCellCallsOutChannelDashboardPlots
    
    script:
    """
    outSummaryPlot="${inCountReportsMap['sample']}-called_cells_summary.pdf"
    outSummaryStats="${inCountReportsMap['sample']}-called_cells_summary.stats.txt"
    
    barnyardParams=""
    if [ "${inBarnyardMap['isBarnyard']}" == 1 ]
    then
        barnyardParams="--window_matrices ${inWindowMatrix} --barnyard"
    fi
    
    Rscript ${script_dir}/summarize_cell_calls.R \
        --sample_name ${inCountReportsMap['sample']} \
        --read_count_tables ${inCountReports} \
        --stats_files ${inCalledCellStats} \
        --insert_size_tables ${inInsertSizes} \
        --peak_call_files ${inNarrowPeaks} \
        --merged_peaks ${inMergedPeaks} \
        --per_base_tss_region_coverage_files ${inPerBaseCoverageTss} \
        --plot \${outSummaryPlot} \
        --output_stats \${outSummaryStats} \${barnyardParams}
    """
}

summarizeCellCallsOutChannelCallCellsSummaryStats
    .into { summarizeCellCallsOutChannelCallCellsSummaryStatsCopy01;
            summarizeCellCallsOutChannelCallCellsSummaryStatsCopy02 }


/*
** ================================================================================
** Make genome browser files.
** ================================================================================
*/

sortTssBedOutChannelCopy03
    .toList()
    .flatMap{ makeGenomeBrowserFilesChannelSetupTss( it, sampleSortedNames, sampleGenomeMap ) }
    .set { makeGenomeBrowserFilesInChannelTssRegions }

mergePeaksByFileOutChannelCopy05
    .toList()
    .flatMap{ makeGenomeBrowserFilesChannelSetupMergedPeaks( it, sampleSortedNames ) }
    .set { makeGenomeBrowserFilesInChannelMergedPeaks }

getUniqueFragmentsOutChannelTranspositionSitesCopy08
    .flatten()
    .toList()
    .flatMap{ makeGenomeBrowserFilesChannelSetupTranspositionSites( it, sampleSortedNames ) }
    .set { makeGenomeBrowserFilesInChannelTranspositionSites }

sortChromosomeSizeOutChannelCopy09
    .toList()
    .flatMap { makeGenomeBrowserFilesChannelSetupChromosomeSizes( it, sampleSortedNames, sampleGenomeMap ) }
    .set { makeGenomeBrowserFilesInChannelChromosomeSizes }

process makeGenomeBrowserFilesProcess {
    cache 'lenient'
    errorStrategy onError
    module 'openjdk/latest:modules:modules-init:modules-gs'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "genome_browser" ) }, pattern: "*.tss_file.sorted.gb.*", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "genome_browser" ) }, pattern: "*-merged_peaks.gb.*", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "genome_browser" ) }, pattern: "*-transposition_sites.gb.*", mode: 'copy'

    input:
    set file( inTssRegions ), inTssRegionsMap from makeGenomeBrowserFilesInChannelTssRegions
    set file( inMergedPeaks ), inMergedPeaksMap from makeGenomeBrowserFilesInChannelMergedPeaks
    set file( inTranspositionSites ), file( inTbiTranspositionSites ), inTranspositionSitesMap from makeGenomeBrowserFilesInChannelTranspositionSites
    set file( inChromosomeSizes ), inChromosomeSizesMap from makeGenomeBrowserFilesInChannelChromosomeSizes

    output:
    file( "*.tss_file.sorted.gb.bed*" ) into makeGenomeBrowserFilesProcessOutChannelTssRegionsBed
    file( "*-merged_peaks.gb.bed" ) into makeGenomeBrowserFilesProcessOutChannelMergedPeaksBed
    file( "*-transposition_sites.gb.bed*" ) into makeGenomeBrowserFilesProcessOutChannelTranspositionSitesBed
    file( "*.tss_file.sorted.gb.bb" ) into makeGenomeBrowserFilesProcessOutChannelTssRegionsBB
    file( "*-merged_peaks.gb.bb" ) into makeGenomeBrowserFilesProcessOutChannelMergedPeaksBB
    file( "*-transposition_sites.gb.bb" ) into makeGenomeBrowserFilesProcessOutChannelTranspositionSitesBB

    script:
    """
    outTssGb="${inTssRegionsMap['sample']}-${inTssRegionsMap['genome']}.tss_file.sorted.gb"
    outMergedPeaksGb="${inMergedPeaksMap['sample']}-merged_peaks.gb"
    outTranspositionSitesGb="${inTranspositionSitesMap['sample']}-transposition_sites.gb"
    
    awk 'BEGIN{OFS="\t"}{\$1="chr" \$1}1' $inChromosomeSizes > chromosome_size.txt.edited

    zcat $inTssRegions | awk 'BEGIN{OFS="\t"}{\$1="chr" \$1}1' | LC_COLLATE=C sort -k 1,1 -k2,2n > \${outTssGb}.bed
    bedToBigBed -tab \${outTssGb}.bed chromosome_size.txt.edited \${outTssGb}.bb
    bgzip \${outTssGb}.bed
    tabix -p bed \${outTssGb}.bed.gz

    cat $inMergedPeaks | awk 'BEGIN{OFS="\t"}{\$1="chr" \$1}1' | LC_COLLATE=C sort -k 1,1 -k2,2n > \${outMergedPeaksGb}.bed
    bedToBigBed -tab \${outMergedPeaksGb}.bed chromosome_size.txt.edited \${outMergedPeaksGb}.bb
    
    zcat $inTranspositionSites | $script_dir/trim_bed_stream.py $inChromosomeSizes | awk 'BEGIN{OFS="\t"}{\$1="chr" \$1}1' | LC_COLLATE=C sort -k 1,1 -k2,2n > \${outTranspositionSitesGb}.bed
    bedToBigBed -tab \${outTranspositionSitesGb}.bed chromosome_size.txt.edited \${outTranspositionSitesGb}.bb
    bgzip \${outTranspositionSitesGb}.bed
    tabix -p bed \${outTranspositionSitesGb}.bed.gz
    """
}


/*
** ================================================================================
** Get per cell insert sizes and banding scores (optional).
** ================================================================================
*/

/*
** Input files.
**   ajh var                 ajh filename                   nf channel
**   fragments_file          %s.fragments.txt.gz            getUniqueFragmentOutChannelFragmentsCopy01
**   cell_whitelist          %s.cell_whitelist.txt          callCellsOutChannelCalledCellsWhitelistCopy04
** Output files.
**   per_cell_insert_sizes   %s.per_cell_insert_sizes.txt   getPerCellInsertSizesOutChannel
**   banding_scores          %s.banding_scores.txt          getBandingScoresOutChannel
*/

getUniqueFragmentOutChannelFragmentsCopy01
	.flatten()
    .toList()
    .flatMap{ getBandingScoresChannelSetupFragments( it, sampleSortedNames ) }
    .set { getBandingScoresInChannelFragments }
    
callCellsOutChannelCalledCellsWhitelistCopy04
    .toList()
    .flatMap { getBandingScoresChannelSetupCellWhitelist( it, sampleSortedNames ) }
    .set { getBandingScoresInChannelCellWhitelist }

process getBandingScoresProcess {
    cache 'lenient'
    errorStrategy onError

    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "banding_scores" ) }, pattern: "*-per_cell_insert_sizes.txt", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "banding_scores" ) }, pattern: "*-banding_scores.txt", mode: 'copy'

    input:
    set file( inFragments ), file(inFragmentsTbi), inFragmentsMap from getBandingScoresInChannelFragments
    set file( inCellWhitelist ), inCellWhiteListMap from getBandingScoresInChannelCellWhitelist
    
    output:
    file( "*-per_cell_insert_sizes.txt" ) into getBandingScoresOutChannelCellInsertSizes
    file( "*-banding_scores.txt" ) into getBandingScoresOutChannelBandingScores

	when:
		params.calculate_banding_scores

    script:
    """
    # output filenames
    source ${pipeline_path}/load_python_env_reqs.sh
    source ${script_dir}/python_env/bin/activate
    
    outPerCellInsertSizesFile="${inFragmentsMap['sample']}-per_cell_insert_sizes.txt"
    outBandingScoresFile="${inFragmentsMap['sample']}-banding_scores.txt"
    
    python ${script_dir}/get_insert_size_distribution_per_cell.py ${inFragments} \${outPerCellInsertSizesFile} --barcodes ${inCellWhitelist}
    
    Rscript ${script_dir}/calculate_nucleosome_banding_scores.R \${outPerCellInsertSizesFile} \${outBandingScoresFile} --barcodes ${inCellWhitelist}
    """
}


/*
** ================================================================================
** Call motifs.
** ================================================================================
*/

/*
** Input files.
**   ajh var                 ajh filename           nf channel
**   fasta  (from genome.json file path)
**   merged_peaks                                   *-merged_peaks.bed
**   motifs (from genome.json file path)
** Output files.
**   output_file/peak_motif_files  gc_binned.%s.peak_calls.bb xxx
** Notes:
**   o  require genome information, which depends on sample
**   o  there is an 'each' input channel thingy, which is likely
**      to accomplish the same function as the .combine() operator
**      below. It may make sense to modify this sometime in the
**      future, in order to simplify this process a bit. See
**      https://www.nextflow.io/docs/latest/faq.html  How do I iterate over a process n times?
**      process bootstrapReplicateTrees {
**        publishDir "$results_path/$datasetID/bootstrapsReplicateTrees"
**
**        input:
**        each x from 1..bootstrapReplicates
**        set val(datasetID), file(ClustalwPhylips)
**
**        output:
**        file "bootstrapTree_${x}.nwk" into bootstrapReplicateTrees
**
**        script:
*/

Channel
	.fromList( 0..params.motif_calling_gc_bins-1 )
	.set { gcBinsInChannel }

mergePeaksByFileOutChannelCopy06
	.combine( gcBinsInChannel )
	.toList()
	.flatMap { callMotifsChannelSetupMergedPeaks( it, sampleSortedNames, sampleGenomeMap ) }
    .set { callMotifsInChannelCombined }

process callMotifsProcess {
	cache 'lenient'
    errorStrategy onError

    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "peak_motifs" ) }, pattern: "*-peak_calls.bb", mode: 'copy'

	input:
	set file( inMergedPeaks ), inMergedPeaksMap from callMotifsInChannelCombined

	output:
	file( "*-peak_calls.bb" ) into callMotifsOutChannelPeakCalls
	
	when:
		inMergedPeaksMap['fasta'] != null && inMergedPeaksMap['motifs'] != null
		
	script:
	"""
	source ${pipeline_path}/load_python_env_reqs.sh
	source ${script_dir}/python_env/bin/activate
	
	gc_bin_padded=`echo ${inMergedPeaksMap['gc_bin']} | awk '{printf("%02d",\$1+1)}'`
	outGcBinned="${inMergedPeaksMap['sample']}-gc_\${gc_bin_padded}-peak_calls.bb"
	
	python ${script_dir}/call_peak_motifs.py ${inMergedPeaksMap['fasta']} ${inMergedPeaks} ${inMergedPeaksMap['motifs']} \${outGcBinned} --gc_bin ${inMergedPeaksMap['gc_bin']} --pwm_threshold ${task.ext.pwm_threshold}
	"""
}


/*
** ================================================================================
** Make motif matrix
** ================================================================================
*/

/*
** Input files.
**   ajh var                 ajh filename           nf channel
**   peak_motif_files        'gc_binned.%s.peak_calls.bb'  ${inMergedPeaksMap['sample']}-gc_\${gc_bin_padded}-peak_calls.bb    (as string of space-separated filenames)
**   fasta (GENOME_FILES)
**   merged_peaks            'merged_peaks.bed'     *-merged_peaks.bed
**   motifs (GENOME_FILES)
** Output files.
**   peak_tf_matrix          'peak_motif_matrix.mtx.gz'  *-peak_motif_matrix.mtx.gz
*/

callMotifsOutChannelPeakCalls
	.toList()
	.flatMap { makeMotifMatrixChannelSetupPeakCalls( it, sampleSortedNames, sampleGenomeMap ) }
	.set { makeMotifMatrixInChannelPeakCalls }
	
mergePeaksByFileOutChannelCopy07
	.toList()
	.flatMap { makeMotifMatrixChannelSetupMergedPeaks( it, sampleSortedNames, sampleGenomeMap ) }
	.set { makeMotifMatrixInChannelMergedPeaks }
	
process makeMotifMatrixProcess {
	cache 'lenient'
    errorStrategy onError

    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "motif_matrices" ) }, pattern: "*-peak_motif_matrix.*", mode: 'copy'

	input:
		tuple file( inPeakCalls ), inPeakCallsMap from makeMotifMatrixInChannelPeakCalls
		tuple file( inMergedPeaks ), inMergedPeaksMap from makeMotifMatrixInChannelMergedPeaks

	output:
	file( "*-peak_motif_matrix.*" ) into makeMotifMatrixOutChannel

	when:
		inPeakCallsMap['fasta'] != null && inPeakCallsMap['motifs'] != null

	script:
	"""
	source ${pipeline_path}/load_python_env_reqs.sh
	source ${script_dir}/python_env/bin/activate
	
	outPeakTfMatrix="${inPeakCallsMap['sample']}-peak_motif_matrix.mtx.gz"
	
	python ${script_dir}/generate_motif_matrix.py \
	--peak_motif_files ${inPeakCalls} \
	--fasta ${inPeakCallsMap['fasta']} \
	--peaks ${inMergedPeaks} \
	--motifs ${inPeakCallsMap['motifs']} \
	--peak_tf_matrix \${outPeakTfMatrix}
	"""
}


/*
** ================================================================================
** Make reduced dimension matrix.
** ================================================================================
*/

/*
** Input files.
**   ajh var                   ajh filename           nf channel
**   peak_matrices             "*-peak_matrix.mtx.gz"
**   promoter_matrices         "*-promoter_matrix.mtx.gz"
** Output files.
**   (tfidf_matrix)            peak_motif_matrix.columns.txt  peak_motif_matrix.mtx.gz  peak_motif_matrix.rows.txt
**   (svd_coords)              '%s.svd_coords.txt'           replace with pca_coords.txt
**   (umap_coords)             '%s.umap_coords.txt'
**   (tsne_coords)             '%s.tsne_coords.txt'
**   (umap_plot)                                             <sample>.umap_plot.pdf
**   (monocle3 object)  '%s.monocle3_cds.rds'
** Notes:
**   o  the reduce_dimensions.R script exits without error when
**      the filtered peak/promoter matrices have no cells. In this
**      case, the script does not create output files. In order for
**      Nextflow to continue despite such errors, set process
**      directive
**
**         "errorStrategy 'ignore'"
**   o  reduce_dimensions.R should exit with 0 status
**   o  scrublet may exit with non-zero status so errorStrategy is 'ignore'
**   o  create empty files in case the programs fail to create the files
**      so that the process block publishes log files
*/

makePeakMatrixOutChannel
	.toList()
	.flatMap { makeReducedDimensionMatrixChannelSetupPeakMatrix( it, sampleSortedNames, sampleGenomeMap ) }
	.set { makeReducedDimensionMatrixInChannelPeakMatrix }
	
makeCountReportsOutChannelCopy03
    .toList()
    .flatMap { makeReducedDimensionMatrixChannelSetupCountReport( it, sampleSortedNames ) }
    .set { makeReducedDimensionMatrixInChannelCountReport }

process makeReducedDimensionMatrixProcess {
	cache 'lenient'
    errorStrategy 'ignore'

    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "reduce_dimension" ) }, pattern: "*-scrublet_table.csv", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "reduce_dimension" ) }, pattern: "*-lsi_coords.txt", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "reduce_dimension" ) }, pattern: "*-umap_coords.txt", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "reduce_dimension" ) }, pattern: "*-umap_plot.*", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "reduce_dimension" ) }, pattern: "*-monocle3_cds.rds", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "reduce_dimension" ) }, pattern: "*-blacklist_regions_file.log", mode: 'copy'
    publishDir path: "${analyze_dir}", saveAs: { qualifyFilename( it, "reduce_dimension" ) }, pattern: "*-reduce_dimensions.log", mode: 'copy'
    publishDir path: "${output_dir}/analyze_dash/img", pattern: "*-umap_plot.png", mode: 'copy'
    publishDir path: "${output_dir}/analyze_dash/img", pattern: "*-scrublet_hist.png", mode: 'copy'

	input:
	tuple file( inPeakFiles ), inPeakMatrixMap from makeReducedDimensionMatrixInChannelPeakMatrix
    tuple file( inCountReport ), inCountReportMap from makeReducedDimensionMatrixInChannelCountReport

	output:
    file("*-lsi_coords.txt") into makeReducedDimensionMatrixOutChannelPcaCoords
    file("*-umap_coords.txt") into makeReducedDimensionMatrixOutChannelUmapCoords
    file("*-umap_plot.*") into makeReducedDimensionMatrixOutChannelUmapPlot
    file("*-monocle3_cds.rds") into makeReducedDimensionMatrixOutChannelMonocle3Cds
    file("*-scrublet_hist.png") into makeReducedDimensionMatrixOutChannelScrubletHist
    file("*-scrublet_table.csv") into makeReducedDimensionMatrixOutChannelScrubletTable
    file("*-blacklist_regions_file.log") into makeReducedDimensionMatrixOutChannelBlackListRegionsFile
    file("*-reduce_dimensions.log") into makeReducedDimensionMatrixOutChannelReducedDimensionsLogFile

	script:
	"""
    inPeakMatrix="${inPeakMatrixMap['sample']}-peak_matrix.mtx.gz"
    inSampleName="${inPeakMatrixMap['sample']}"

   	 outScrubletHistFile="${inPeakMatrixMap['sample']}-scrublet_hist.png"
    outScrubletTableFile="${inPeakMatrixMap['sample']}-scrublet_table.csv"
	outLsiCoordsFile="${inPeakMatrixMap['sample']}-lsi_coords.txt"
	outUmapCoordsFile="${inPeakMatrixMap['sample']}-umap_coords.txt"
	outUmapPlotFile="${inPeakMatrixMap['sample']}-umap_plot"
	outMonocle3CdsFile="${inPeakMatrixMap['sample']}-monocle3_cds.rds"
    outBlackListRegionsFile="${inPeakMatrixMap['sample']}-blacklist_regions_file.log"

    umi_cutoff=$task.ext.umi_cutoff
    frip_cutoff=$task.ext.frip_cutoff
    frit_cutoff=$task.ext.frit_cutoff
    num_lsi_dimensions=$task.ext.num_lsi_dimensions
    cluster_resolution=$task.ext.cluster_resolution
	doublet_predict_top_ntile=$task.ext.doublet_predict_top_ntile

    doublet_predict=""
    if [ "${params.doublet_predict}" ]
    then
        doublet_predict=" --doublet_predict "
    fi

    black_list_file=""
    echo "blacklist_regions: ${inPeakMatrixMap['blacklist_regions']}"

    if [[ "${params.filter_blacklist_regions}" == "true" && "${inPeakMatrixMap['blacklist_regions']}" != "" ]]
    then
        black_list_file=" --black_list_file ${inPeakMatrixMap['blacklist_regions']} "
        echo "blacklist_region file: ${inPeakMatrixMap['blacklist_regions']}" > \${outBlackListRegionsFile}
    else
      echo "blacklist_region file: none" > \${outBlackListRegionsFile}
    fi

    if [ "${params.doublet_predict}" ]
    then
        $script_dir/run_scrublet.py --sample_name=\${inSampleName} --mat_file=\${inPeakMatrix} --umi_cutoff=\$umi_cutoff
    fi

	Rscript ${script_dir}/reduce_dimensions.R \
    --sample_name \${inSampleName} \
	--mat_file \${inPeakMatrix} \
    --count_file ${inCountReport} \
    --umi_cutoff \${umi_cutoff} \
    --frip_cutoff \${frip_cutoff} \
    --frit_cutoff \${frit_cutoff} \
    --doublet_predict_top_ntile \${doublet_predict_top_ntile} \
    --num_lsi_dimensions \${num_lsi_dimensions} \
    --cluster_resolution \${cluster_resolution} \
    --cds_file \${outMonocle3CdsFile} \
    --lsi_coords_file \${outLsiCoordsFile} \
    --umap_coords_file \${outUmapCoordsFile} \
    --umap_plot_file \${outUmapPlotFile} \${doublet_predict} \${black_list_file}

    if [ ! -e "\${outScrubletHistFile}" ]
    then
      touch "\${outScrubletHistFile}"
    fi

    if [ ! -e "\${outScrubletTableFile}" ]
    then
      touch "\${outScrubletTableFile}"
    fi

    if [ ! -e "\${outLsiCoordsFile}" ]
    then
      touch "\${outLsiCoordsFile}"
    fi

    if [ ! -e "\${outUmapCoordsFile}" ]
    then
      touch "\${outUmapCoordsFile}"
    fi

    if [ ! -e "\${outUmapPlotFile}.png" ]
    then
      touch "\${outUmapPlotFile}.null"
    fi

    if [ ! -e "\${outMonocle3CdsFile}" ]
    then
      touch "\${outMonocle3CdsFile}"
    fi

    if [ ! -e "\${outBlackListRegionsFile}" ]
    then
      touch "\${outBlackListRegionsFile}"
    fi
	"""
}


/*
** ================================================================================
** Make experiment dashboard.
** ================================================================================
*/

summarizeCellCallsOutChannelCallCellsSummaryStatsCopy01
	.toList()
	.map { experimentDashboardProcessChannelSetup( it, sampleSortedNames, sampleGenomeMap ) }
	.set { experimentDashboardProcessInChannel }
	
process experimentDashboardProcess {
	cache 'lenient'
    errorStrategy onError
	publishDir path: "${output_dir}/analyze_dash/js", pattern: "run_data.js", mode: 'copy'
	
	input:
	file( "*") from experimentDashboardProcessInChannel

	output:
	file( "run_data.js" ) into experimentDashboardProcessOutChannelRunData

	script:
	"""
	mkdir -p ${output_dir}/analyze_dash/js
	${script_dir}/make_run_data.py -i ${output_dir}/analyze_out/args.json -o run_data.js
	
	mkdir -p ${output_dir}/analyze_dash/js ${output_dir}/analyze_dash/img
	cp ${script_dir}/skeleton_dash/img/*.png ${output_dir}/analyze_dash/img
	cp ${script_dir}/skeleton_dash/js/* ${output_dir}/analyze_dash/js
	cp -r ${script_dir}/skeleton_dash/style ${output_dir}/analyze_dash
	cp ${script_dir}/skeleton_dash/exp_dash.html ${output_dir}/analyze_dash
	"""
}


/*
** ================================================================================
** Make merged (sample) plot files.
** ================================================================================
*/

summarizeCellCallsOutChannelCallCellsSummaryStatsCopy02
    .toList()
    .map { makeMergedPlotFilesProcessChannelSetupCallCellsSummaryStats( it, sampleSortedNames, sampleGenomeMap ) }
    .set { makeMergedPlotFilesProcessInChannelCallCellsSummaryStats }

makeReducedDimensionMatrixOutChannelUmapPlot
    .toList()
    .map { makeMergedPlotFilesProcessChannelSetupMakeMergedUmapPlots( it, sampleSortedNames, sampleGenomeMap ) }
    .set { makeMergedPlotFilesProcessInChannelMakeMergedUmapPlots }

process makeMergedPlotFilesProcess {
    cache 'lenient'
    errorStrategy onError
    publishDir path: "${output_dir}/analyze_out/merged_plots", pattern: "merged.called_cells_summary.stats.tsv", mode: 'copy'
    publishDir path: "${output_dir}/analyze_out/merged_plots", pattern: "merged.called_cells_summary.pdf", mode: 'copy'
    publishDir path: "${output_dir}/analyze_out/merged_plots", pattern: "merged.umap_plots.pdf", mode: 'copy'

    input:
    file( "*") from makeMergedPlotFilesProcessInChannelCallCellsSummaryStats
    file( "*-umap_plots.pdf") from makeMergedPlotFilesProcessInChannelMakeMergedUmapPlots

    output:
    file( "merged.called_cells_summary.stats.tsv" ) into makeMergedPlotFilesProcessOutChannelMergedCalledCellsSummaryTsv
    file( "merged.called_cells_summary.pdf" ) into makeMergedPlotFilesProcessOutChannelMergedCalledCellsSummaryPdf
    file( "merged.umap_plots.pdf" ) into makeMergedPlotFilesProcessOutChannelMergedUmapPlots

    script:
    """
    mkdir -p ${output_dir}/analyze_out/merged_plots
    ${script_dir}/merge_summary_plots.py -i ${output_dir}/analyze_out/args.json -o merged.called_cells_summary.pdf
    ${script_dir}/merge_umap_plots.py -i ${output_dir}/analyze_out/args.json -o merged.umap_plots.pdf

    header='sample cell_threshold fraction_hs fraction_tss median_per_cell_frip median_per_cell_frit tss_enrichment sample_peaks_called total_merged_peaks total_reads fraction_reads_in_cells total_barcodes number_of_cells median_reads_per_cell min_reads_per_cell max_reads_per_cell median_duplication_rate median_fraction_molecules_observed median_total_fragments total_deduplicated_reads [bloom_collision_rate]'
    stats_file='merged.called_cells_summary.stats.tsv'
    header_wtabs=`echo \${header} | sed 's/ /\t/g'`
    rm -f \${stats_file}
    echo "\${header_wtabs}" > \${stats_file}
    lfil=`ls *-called_cells_summary.stats.txt`
    for fil in \${lfil}
    do
      tail -n +2 \${fil} >> \${stats_file}
    done
    """
}

/*
** Template
**
process template {
	cache 'lenient'
    errorStrategy onError

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


def printErr( errString ) {
	System.err.println( errString )
}


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
    log.info '    --help                                           Show this message and exit.'
    log.info ''
    log.info 'Required parameters (specify in your config file):'
    log.info '    params.run_dir = RUN_DIRECTORY (str)             Path to the Illumina run directory.'
    log.info '    params.sample_sheet = SAMPLE_SHEET (str)         Path to sample sheet file.'
    log.info '    params.output_dir = OUTPUT_DIRECTORY (str)       Path to the demux directory with the demultiplexed fastq files.'
    log.info '    params.genomes_json = GENOMES_JSON               A json file of genome information for analyses.'
    log.info ''
    log.info 'Optional parameters (specify in your config file):'
    log.info '    params.bowtie_cpus = CPUS (int)                  Number of threads in bowtie run.'
    log.info '    params.bowtie_memory = MEMORY (int)              Amount of memory in bowtie run.'
    log.info '    params.bowtie_seed = SEED (int)                  Bowtie random number generator seed.'
    log.info '    params.doublet_predict VALUE (logical)           Run double prediction.'
    log.info '    params.filter_blacklist_regions VALUE (logical)  Filter blacklisted genomic regions when preprocessing the cell_data_set.'
	log.info '    params.calculate_banding_scores (flag)           Add flag to calculate banding scores (take a while to run which is why opt-in here).'
    log.info '    params.samples = SAMPLES (str)                   A list of samples to include in the analysis. Restrict to a subset of samples in demux output.'
    log.info '    params.reads_threshold = THRESHOLD (int)         Threshold for how many unique reads to consider a cell. Default is to dynamically pick with mclust.'
    log.info '    process.max_forks = 20                           The maximum number of processes to run at the same time on the cluster.'
    log.info '    process.queue = "trapnell-short.q"               The queue on the cluster where the jobs should be submitted.'
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
	s += String.format( "Demultiplexed fastq directory:        %s\n", demux_dir )
	s += String.format( "Analysis output directory:            %s\n", analyze_dir )
    s += String.format( "Launch directory:                     %s\n", workflow.launchDir )
    s += String.format( "Work directory:                       %s\n", workflow.workDir )
	if( params.samples != null ) {
		s += String.format( "Samples to include in analysis:   %s\n", params.samples )
	}
    if( params.bowtie_seed != 0 ) {
        s += String.format( "Bowtie random number generator seed:  %s\n", params.bowtie_seed )
    }
	if( params.reads_threshold != null ) {
	   s += String.format( "Threshold for # reads/cell:           %s\n", params.reads_threshold )
	}
	if( params.calculate_banding_scores != null ) {
	   s += String.format( "Calculate banding scores:             %s\n", params.calculate_banding_scores )
	}
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
** Function: writeChannelEntries
** Purpose: diagnostic
** Use: place it within a channel transformation 'chain' using the .map operator.
** For example,
**
**     runAlignOutChannel
**	  .toList()
**	  .map { writeChannelEntries( it, 'runAlignOutChannel start' ) }
**	  .flatMap { mergeBamChannelSetup( it, sampleLaneMap ) }
**	  .map { writeChannelEntries( it, 'runAlignOutChannel finish' ) }
**    .set { mergeBamsInChannel }
*/
def writeChannelEntries( inEntry, label ) {
	println( "Channel entries: $label: $inEntry" )
	println()
	return( inEntry )
}


/*
** Make directory, if it does not exist.
*/
def makeDirectory( directoryName ) {
	dirHandle = new File( directoryName )
	if( !dirHandle.exists() ) {
		if( !dirHandle.mkdirs() ) {
			printErr( "Error: unable to create directory $directoryName" )
			System.exit( -1 )
		}
	}
}


/*
** Check that directory exists and can be read.
*/
def checkDirectory( directoryName ) {
	def dirHandle = new File( directoryName )
	if( !dirHandle.exists() ) {
		printErr( "Error: unable to find directory $directoryName" )
		return( false )
	}
	if( !dirHandle.canRead() ) {
		printErr( "Error: unable to read directory $directoryName" )
		return( false )
	}
	return( true )
}


/*
** Check that file exists and can be read.
*/
def checkFile( fileName ) {
	def fileHandle = new File( fileName )
	if( !fileHandle.exists() ) {
		printErr( "Error: unable to find file $fileName" )
		return( false )
	}
	if( !fileHandle.canRead() ) {
		printErr( "Error: unable to read file $fileName" )
		return( false )
	}
	
	return( true )
}


/*
** Check directories for existence and/or accessibility.
*/
def checkDirectories( params ) {
	def dirName
	
	/*
	** Check that demux_dir exists.
	*/
	dirName = demux_dir
	if( !checkDirectory( dirName ) ) {
		System.exit( -1 )
	}
	
	/*
	** Check that either the analyze_dir exists or we can create it.
	*/
	dirName = analyze_dir
	makeDirectory( dirName )
}


/*
** Make archive copies of various files and store in reports sub-directory.
*/
def archiveRunFiles( params, timeNow )
{
	makeDirectory( "${analyze_dir}/reports" )
	file_suffix = timeNow.format( 'yyyy-MM-dd_HH-mm-ss' )
	def i = 1
	Path dst
	workflow.configFiles.each { aFile ->
    src = aFile
    dst = Paths.get( "${analyze_dir}/reports/${aFile.getName()}.${file_suffix}.${i}" )
    Files.copy( src, dst )
    i += 1
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
        if( samples == null ) {
            printErr "Error: no samples list for run " + aRun + " in args.json"
            System.exit( -1 )
        }
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


def getSamplePeakGroupMap( argsJson ) {
    /*
    ** Make a map of peak groups keyed by sample. And check
    ** that the groups are consistent, if there is more than
    ** one run. Allow for empty, zero-length, values.
    */
    def peakGroups
    def samplePeakGroupMap = [:]
    def errorFlag = 0
    def runs = argsJson.keySet()
    runs.each { aRun ->
        peakGroups = argsJson[aRun]['peak_groups']
        if( peakGroups == null ) {
            printErr "Error: no peak_groups map for run " + aRun + " in args.json"
            System.exit( -1 )
        }
        peakGroups.each { aSample, aGroup ->
            if( ! samplePeakGroupMap.containsKey( aSample ) ) {
                samplePeakGroupMap[aSample] = aGroup
            } else {
                if( samplePeakGroupMap[aSample] != aGroup ) {
                    printErr "Error: inconsistent peak groups assigned to sample: " + aSample
                    errorFlag = 1
                }
            }
        }
    }
    if( errorFlag == 1 ) {
        System.exit( -1 )
    }
    return( samplePeakGroupMap )
}


def getSamplePeakFileMap( argsJson ) {
    /*
    ** Make a map of peak files keyed by sample. And check
    ** that the files are consistent, if there is more than
    ** one run. Allow for empty, zero-length, values.
    */
    def peakFiles
    def samplePeakFileMap = [:]
    def errorFlag = 0
    def runs = argsJson.keySet()
    runs.each { aRun ->
        peakFiles = argsJson[aRun]['peak_files']
        if( peakFiles == null ) {
            printErr "Error: no peak_files map for run " + aRun + " in args.json"
            System.exit( -1 )
        }
        peakFiles.each { aSample, aFile ->
            if( ! samplePeakFileMap.containsKey( aSample ) ) {
                samplePeakFileMap[aSample] = aFile
            } else {
                if( samplePeakFileMap[aSample] != aFile ) {
                    printErr "Error: inconsistent peak files assigned to sample: " + aSample
                    errorFlag = 1
                }
            }
        }
    }
    if( errorFlag == 1 ) {
        System.exit( -1 )
    }
    return( samplePeakFileMap )
}


/*
** Check that sample-specific peak files exist.
*/
def checkPeakFiles( samplePeakFileMap ) {
    def errorFlag = 0
    def fileHandle
    samplePeakFileMap.each { aSample, aFile ->
        if( aFile.length() > 0 ) {
            fileHandle = new File( aFile )
            if( !fileHandle.exists() ) {
                printErr "Error: unable to find peak file \'${aFile}\'"
                errorFlag = 1
            }
            if( !fileHandle.canRead() ) {
                printErr "Error: unable to read peak file \'${aFile}\'"
                errorFlag = 1
            }
        }
    }
    if( errorFlag == 1 ) {
        System.exit( -1 )
    }
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
	def dirName = demux_dir
	def samples = sampleLaneMap.keySet()
	samples.each { aSample ->
		def lanes = sampleLaneMap[aSample]
		lanes.each { aLane ->
			fileName = dirName + '/' + aSample + '/fastqs_trim/' + String.format( '%s-%s_R1.trimmed.fastq.gz', aSample, aLane )
			fileHandle = new File( fileName )
			if( !fileHandle.exists() ) {
				printErr( "Error: unable to find fastq file ${fileName}" )
				System.exit( -1 )
			}
			if( !fileHandle.canRead() ) {
				printErr( "Error: unable to read fastq file ${fileName}" )
				System.exit( -1 )
			}
			fileName = dirName + '/' + aSample + '/fastqs_trim/' + String.format( '%s-%s_R2.trimmed.fastq.gz', aSample, aLane )
			fileHandle = new File( fileName )
			if( !fileHandle.exists() ) {
				printErr( "Error: unable to find fastq file ${fileName}" )
				System.exit( -1 )
			}
			if( !fileHandle.canRead() ) {
				printErr( "Error: unable to read fastq file ${fileName}" )
				System.exit( -1 )
			}
		}
	}
}


/*
** Get sample names from demux args.json file and return in a list.
*/
def getArgsJsonSamples( args_json ) {
	def samples = []
	args_json.each { key, run ->
		run['samples'].each { sample ->
			samples.add( sample )
		}
	}
	return( samples.unique() )
}


/*
** Write run data JSON file(s).
*/
def writeRunDataJsonFile( params, argsJson, sampleGenomeMap, jsonFilename, timeNow ) {
	def samples
	if( params.samples ) {
		samples = params.samples
	}
	else {
		samples = getArgsJsonSamples( argsJson )
	}
	
    analyzeDict = [:]
    analyzeDict['run_date'] = timeNow.format( 'yyyy-MM-dd_HH-mm-ss' )
    analyzeDict['demux_dir'] = demux_dir
    analyzeDict['analyze_dir'] = analyze_dir
    analyzeDict['genomes_json'] = params.genomes_json
    analyzeDict['samples'] = samples
    analyzeDict['bowtie_seed'] = params.bowtie_seed
    analyzeDict['reads_threshold'] = params.reads_threshold
    analyzeDict['calculate_banding_scores'] = params.calculate_banding_scores
    analyzeDict['RUNS'] = argsJson
    File file_json = new File( jsonFilename )
    file_json.write( JsonOutput.prettyPrint( JsonOutput.toJson( analyzeDict ) ) )
}


/*
** Check for required genomes.
** Notes:
**   o  perhaps check that the required genome files exist...
*/
def readGenomesJson( params, argsJson ) {
	def pathFil = params.genomes_json
	def fileHandle = new File( pathFil )
	if( !fileHandle.exists() ) {
		printErr( "Error: unable to find genomes json file \'${pathFil}\'" )
		System.exit( -1 )
	}
	if( !fileHandle.canRead() ) {
		printErr( "Error: unable to read genomes json file \'${pathFil}\'" )
		System.exit( -1 )
	}
	def jsonSlurper  = new JsonSlurper()
	def genomesJson = jsonSlurper.parse( fileHandle )
	def namGenomesJson = genomesJson.keySet()
	def runs = argsJson.keySet()
	runs.each { aRun ->
		def argsGenomes = argsJson[aRun]['genomes'].values()
		argsGenomes.each { aGenome ->
			if( !fileHandle.exists() ) {
				printErr( "Error: unable to find genomes json file \'${pathFil}\'" )
				System.exit( -1 )
			}
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
        } else {
            inBed = ''
        }
        tssBedMaps.add( [ 'sample': aSample, 'genome': aGenomeJsonMap['name'], 'hasTss': hasTss, 'inBed': inBed ] )	   
	}
	
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
        } else {
            inTxt = ''
        }
        chromosomeSizeMaps.add( [ 'sample': aSample, 'genome': aGenomeJsonMap['name'], 'hasChromosomeSizes': hasChromosomeSizes, 'inTxt': inTxt ] )
    }

	return( chromosomeSizeMaps )
}


/*
** Set up peak files specified in samplePeakFileMap.
*/
def peakFileChannelSetup( sampleSortedNames, samplePeakFileMap ) {
    def peakFileMap = []
    def fileName
    def f
    sampleSortedNames.each { aSample ->
        if( samplePeakFileMap[aSample].length() > 0 ) {
            f = new File( samplePeakFileMap[aSample] )
            nameBed = f.getName()
            peakFileMap.add( [ 'sample': aSample, 'inBed': samplePeakFileMap[aSample], 'nameBed': nameBed ] )
        }
    }
    return( peakFileMap )
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
			if( !samplesJson.contains( aSample ) ) {
				printErr( "Error: sample \'${aSample}\' in params.samples is not in the args.json file" )
				System.exit( -1 )
			}
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
	def demuxDir = demux_dir
	def seed = ''
	if( params.bowtie_seed != null ) {
	   seed = "--seed ${params.bowtie_seed}"
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
			alignMaps.add( [ 'sample': aSample, 'lane': aLane, 'fastq1':fastq1, 'fastq2':fastq2, 'genome_index':genome_index, 'whitelist':whitelist, 'seed': seed ] )
		}
	}
	
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected  file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
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
		def tuple = new Tuple( sampleLaneBamMap[aSample], [ 'sample': aSample ] )
		outTuples.add( tuple )
	}
	
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
def getUniqueFragmentsChannelSetupBam( inPaths, sampleSortedNames ) {
	/*
	** Check for expected BAM files.
	*/
	def filesExpected = []
	sampleSortedNames.each { aSample ->
		def fileName = aSample + '-merged.bam'
		filesExpected.add( fileName )
	}
    def filesFound = []
    inPaths.each { aPath ->
        filesFound.add( aPath.getFileName().toString() )
    }
    filesExpected.each { aFile ->
    	if( !( aFile in filesFound ) ) {
    		printErr( "Error: missing expected file \'${aFile}\' in channel" )
    		System.exit( -1 )
    	}
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
    inPaths.each { aPath ->
        def aFile = aPath.getFileName().toString()
        def aSample = aFile.split( '-' )[0]
        if( aFile =~ /merged[.]bam$/ ) {
            pathMap[aSample]['bam'] = aPath
        } else if( aFile =~ /merged[.]bam[.]bai$/ ) {
            pathMap[aSample]['bai'] = aPath
        } else {
            println "Warning: getUniqueFragmentsChannelSetupBam: unexpected file \'${fileName}\'"
        }
    }
    
	/*
	** Set up output channel tuples.
	*/
    
    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inBam = pathMap[aSample]['bam']
        def inBai = pathMap[aSample]['bai']
        def tuple = new Tuple( inBam, inBai, [ 'sample':aSample ] )
        outTuples.add( tuple )
    }
        
	return( outTuples )
}


/*
** Set up channel for making read alignment bigwig file.
*/
def getUniqueFragmentsChannelSetupChromosomeSizes( inPaths, sampleSortedNames, sampleGenomeMap ) {
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
	}

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def aGenomeJsonMap = sampleGenomeMap[aSample]
        def inGenomicIntervals = aSample + '-' + aGenomeJsonMap['name'] + '.chromosome_sizes.sorted.txt'
        def tuple = new Tuple( fileMap[inGenomicIntervals], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

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
	inPaths.each { aPath ->
		filesFound.add( aPath.getFileName().toString() )
	}
	filesExpected.each { aFile ->
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}    
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
	inPaths.each { aPath ->
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
	
	def callPeaksTuples = []
	samples.each { aSample ->
		def macs_genome = sampleGenomeMap[aSample]['macs_genome']
		def tuple = new Tuple( pathMap[aSample]['bed'], pathMap[aSample]['tbi'], [ 'sample': aSample, 'macs_genome': macs_genome ] )
		callPeaksTuples.add( tuple )
	}
	
	return( callPeaksTuples )
}


/*
** ================================================================================
** Make merged peak file channel setup functions.
** ================================================================================
*/

/*
** Set up channel of merged peak by group files for downstream
** analysis. Return a list of tuples, a tuple for each peak group.
** Notes:
**   o  this function takes sample-based paths from the input
**      channel
**   o  this functino returns group-based paths to the output
**      channel
**   o  this function identifies samples with groups and 
**      returns maps that map groups to samples and bed
**      file names
*/
def makePeakByGroupFileChannelSetup( inPaths, sampleSortedNames, samplePeakGroupMap ) {
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}    
	}

    /*
    ** Make a map of file path lists keyed by group. This
    ** forms the input channel to the process. Allow for
    ** empty, zero-length, values.
    */
    def groupPaths = [:]
    inPaths.each { aPath ->
        def fileName = aPath.getFileName().toString()
        def aSample = fileName.split( "-" )[0]
        /*
        ** Skip if sample is not assigned to a
        ** peak group.
        */
        if( samplePeakGroupMap[aSample].length() > 0 ) {
            if( !groupPaths.containsKey( samplePeakGroupMap[aSample] ) ) {
                groupPaths[samplePeakGroupMap[aSample]] = []
            }
            groupPaths[samplePeakGroupMap[aSample]].add( aPath )
        }
    }

    /*
    ** Make a map of lists of sample names keyed by group.
    ** Allow for empty, zero length, strings.
    */
    def outSampleLists = [:]
    samplePeakGroupMap.each { aSample, aGroup ->
        /*
        ** Skip if sample is not assigned to a
        ** peak group.
        */
        if( aGroup.length() > 0 ) {
            if( ! outSampleLists.containsKey( aGroup ) ) {
                outSampleLists[aGroup] = []
            }
            outSampleLists[aGroup].add( aSample )
        }
    }

    /*
    ** Make a list of bed filenames by group.
    ** Notes:
    **   o  all sampleSortedNames has all samples
    **   o  samplePeakGroupMap has samples in peak
    **      groups
    **   o  there can be empty listBeds, which become
    **      empty strings
    */
    def listBeds = [:]
    groupPaths.each { aGroup, aList ->
        listBeds[aGroup] = []
        aList.each { aPath ->
           listBeds[aGroup].add( aPath.getFileName().toString() )
        }
    }

    /*
    ** Make channel list.
    */
    def outTuples = []
    groupPaths.each { aGroup, aList ->
        def tuple = new Tuple( aList, [ 'group' : aGroup,
                               'outSampleList' : outSampleLists[aGroup].join( ' ' ),
                               'listBeds': listBeds[aGroup].join( ' ' ) ] )
        outTuples.add( tuple )
    }

	return( outTuples )
}


/*
** Three possibilities exist:
**   o  a group is given for the sample in the
**      samplesheet file
**   o  an external bed file is given for the sample
**      in the samplesheet file
**   o  both a group and an external bed file are
**      given for the sample
**   o  the inPaths contains the combined group bed
**      and external bed paths so there is at least
**      one path per sample
*/
def makePeakByFileFileChannelSetup( inPaths, sampleSortedNames, samplePeakGroupMap, samplePeakFileMap ) {
    /*
    ** Check for expected files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        if( samplePeakGroupMap[aSample].length() > 0 ) {
            def fileName = aSample + '-group_merged_peaks.bed'
            filesExpected.add( fileName )
        }
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
        if( !( aFile in filesFound ) ) {
            printErr( "Error: : missing expected file \'${aFile}\' in channel" )
            System.exit( -1 )
        }
    }

    /*
    ** Initialize the file paths map. There must
    ** be at least one bed file for each sample.
    */
    filePaths = [:]
    sampleSortedNames.each { aSample ->
        filePaths[aSample] = []
    }

    /*
    ** The merged peak bed files from the merge peak
    ** by group process above.
    */
    inPaths.each { aPath ->
        def fileName = aPath.getFileName().toString()
        def aSample = fileName.split( "-" )[0]
        filePaths[aSample].add( aPath )
    }

    /*
    ** Make a list of bed filenames by sample.
    */
    def listBeds = [:]
    filePaths.each { aSample, aList ->
        listBeds[aSample] = []
        aList.each { aPath ->
            listBeds[aSample].add( aPath.getFileName().toString() )
        }
    }

    /*
    ** Make output channel list.
    */
    def outTuples = []
    filePaths.each { aSample, aList ->
        def tuple = new Tuple( aList, [ 'sample' : aSample,
                                        'listBeds': listBeds[aSample].join( ' ' ) ] )
        outTuples.add( tuple )
    }

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
def makeWindowedGenomeIntervalsChannelSetup( inPaths, sampleSortedNames, sampleGenomeMap ) {
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}    
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def aGenomeJsonMap = sampleGenomeMap[aSample]
        def inGenomicIntervals = aSample + '-' + aGenomeJsonMap['name'] + '.chromosome_sizes.sorted.txt'
        def tuple = new Tuple( fileMap[inGenomicIntervals], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

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
def makePromoterSumIntervalsChannelSetup( inPaths, sampleSortedNames, sampleGenomeMap ) {
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}    
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
        } else {
            hasGeneScoreBed = 0
            inBedName = aSample + '-merged_peaks.bed'
            inBedPath = fileMap[inBedName]
            inBedFile = 'NA'
            tssFile = sampleGenomeMap[aSample]['tss']
        }
        def tuple = new Tuple( inBedPath, [ 'sample': aSample,
                                            'hasGeneScoreBed': hasGeneScoreBed,
                                            'tssFile': tssFile,
                                            'inBedFile': inBedFile ] )
        outTuples.add( tuple )
    }
    
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
    inPaths.each { aPath ->
        filesFound.add( aPath.getFileName().toString() )
    }
    filesExpected.each { aFile ->
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}    
   }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPath ->
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

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def tuple = new Tuple( pathMap[aSample]['bed'], pathMap[aSample]['tbi'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}    
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inMergedPeaks = aSample + '-merged_peaks.bed'
        def tuple = new Tuple( fileMap[inMergedPeaks], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}  
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inFileName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.chromosome_sizes.sorted.txt'
        def tuple = new Tuple( fileMap[inFileName], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

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
    inPaths.each { aPath ->
        filesFound.add( aPath.getFileName().toString() )
    }
    filesExpected.each { aFile ->
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}   
    }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPath ->
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

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def tuple = new Tuple( pathMap[aSample]['bed'], pathMap[aSample]['tbi'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}    
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inBedName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.tss_file.sorted.bed.gz'
        def tuple = new Tuple( fileMap[inBedName], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}    
   }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inFileName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.chromosome_sizes.sorted.txt'
        def tuple = new Tuple( fileMap[inFileName], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }
    
    def outTuples = []
    sampleSortedNames.each { aSample ->
       def inTxt = aSample + '-duplicate_report.txt'
       def tuple = new Tuple( fileMap[inTxt], [ 'sample': aSample ] )
       outTuples.add( tuple )
    }

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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}    
    }
    
    def outTuples = []
    sampleSortedNames.each { aSample ->
       def inTxt = aSample + '-peak_counts.txt'
       def tuple = new Tuple( fileMap[inTxt], [ 'sample': aSample ] )
       outTuples.add( tuple )
    }

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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }
    
    def outTuples = []
    sampleSortedNames.each { aSample ->
       def inTxt = aSample + '-tss_counts.txt'
       def tuple = new Tuple( fileMap[inTxt], [ 'sample': aSample ] )
       outTuples.add( tuple )
    }

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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }
    
    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inTxt = aSample + '-count_report.txt'
        def tuple = new Tuple( fileMap[inTxt], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

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
    inPaths.each { aPath ->
        filesFound.add( aPath.getFileName().toString() )
    }
    filesExpected.each { aFile ->
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPath ->
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

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def tuple = new Tuple( pathMap[aSample]['bed'], pathMap[aSample]['tbi'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inBedName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.tss_file.sorted.bed.gz'
        def tuple = new Tuple( fileMap[inBedName], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inFileName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.chromosome_sizes.sorted.txt'
        def tuple = new Tuple( fileMap[inFileName], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

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
    inPaths.each { aPath ->
        filesFound.add( aPath.getFileName().toString() )
    }
    filesExpected.each { aFile ->
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPath ->
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

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def tuple = new Tuple( pathMap[aSample]['bed'], pathMap[aSample]['tbi'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: unable to find genomes json file \'${aFile}\'" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inMergedPeaks = aSample + '-merged_peaks.bed'
        def tuple = new Tuple( fileMap[inMergedPeaks], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inCellWhitelist = aSample + '-called_cells_whitelist.txt'
        def tuple = new Tuple( fileMap[inCellWhitelist], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
    return( outTuples )
}


/*
** Set up a channel for filtering sites in making peak matrices.
** This channel has the chromosome sizes file.
*/
def makePeakMatrixChannelSetupChromosomeSizes( inPaths, sampleSortedNames, sampleGenomeMap ) {
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inFileName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.chromosome_sizes.sorted.txt'
        def tuple = new Tuple( fileMap[inFileName], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

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
    inPaths.each { aPath ->
        filesFound.add( aPath.getFileName().toString() )
    }
    filesExpected.each { aFile ->
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPath ->
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

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def tuple = new Tuple( pathMap[aSample]['bed'], pathMap[aSample]['tbi'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inGenomicIntervals = aSample + '-genomic_windows.bed'
        def tuple = new Tuple( fileMap[inGenomicIntervals], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inCellWhitelist = aSample + '-called_cells_whitelist.txt'
        def tuple = new Tuple( fileMap[inCellWhitelist], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
    return( outTuples )
}


/*
** Set up a channel for filtering sites in making window matrices.
** This channel has the chromosome sizes file.
*/
def makeWindowMatrixChannelSetupChromosomeSizes( inPaths, sampleSortedNames, sampleGenomeMap ) {
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inFileName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.chromosome_sizes.sorted.txt'
        def tuple = new Tuple( fileMap[inFileName], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

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
    inPaths.each { aPath ->
        filesFound.add( aPath.getFileName().toString() )
    }
    filesExpected.each { aFile ->
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPath ->
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

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def tuple = new Tuple( pathMap[aSample]['bed'], pathMap[aSample]['tbi'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    return( outTuples )
}


def makePromoterMatrixChannelSetupGeneRegions( inPaths, sampleSortedNames, sampleGenomeMap ) {
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inGenomicIntervals = aSample + '-gene_regions.bed.gz'
        def inGeneBodiesGeneMap = sampleGenomeMap[aSample]['gene_bodies_gene_map']
        def inPromoterMatrixRows = aSample + '-promoter_matrix.rows.txt'
        def tuple = new Tuple( fileMap[inGenomicIntervals], [ 'sample': aSample,
                                                              'inGeneBodiesGeneMap': inGeneBodiesGeneMap,
                                                              'inPromoterMatrixRows': inPromoterMatrixRows ] )
        outTuples.add( tuple )
    }
    
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inCellWhitelist = aSample + '-called_cells_whitelist.txt'
        def tuple = new Tuple( fileMap[inCellWhitelist], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
    return( outTuples )
}


/*
** Set up a channel for filtering sites in making promoter matrices.
** This channel has the chromosome sizes file.
*/
def makePromoterMatrixChannelSetupChromosomeSizes( inPaths, sampleSortedNames, sampleGenomeMap ) {
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inFileName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.chromosome_sizes.sorted.txt'
        def tuple = new Tuple( fileMap[inFileName], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inCountReport = aSample + '-count_report.txt'
        def tuple = new Tuple( fileMap[inCountReport], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
    return( outTuples )
}


def summarizeCellCallsSetupCalledCellStats( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected input paths.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-called_cells_stats.json'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inCalledCellStats = aSample + '-called_cells_stats.json'
        def tuple = new Tuple( fileMap[inCalledCellStats], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inInsertSizes = aSample + '-insert_sizes.txt'
        def tuple = new Tuple( fileMap[inInsertSizes], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inCalledPeaks = aSample + '-peaks.narrowPeak.gz'
        def tuple = new Tuple( fileMap[inCalledPeaks], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inMergedPeaks = aSample + '-merged_peaks.bed'
        def tuple = new Tuple( fileMap[inMergedPeaks], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inPerBaseCoverageTss = aSample + '-tss_region_coverage.txt.gz'
        def tuple = new Tuple( fileMap[inPerBaseCoverageTss], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
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
    inPaths.each {  aPath ->
        filesFound.add( aPath.getFileName().toString() )
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPath ->
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

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inWindowMatrix = aSample + '-window_matrix.mtx.gz'
        def tuple = new Tuple( pathMap[aSample]['mtx'], pathMap[aSample]['row'], pathMap[aSample]['col'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }
    
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
    
    return( outMaps )
}


/*
** Set up a channel for getting TSS regions genome browser file.
** This channel has the Tss regions in the genome.
** Notes:
**   o  inPath filenames have the form <sample_name>-<genome_name>.tss_file.sorted.bed.gz
**   o  return a list of files sorted by sample name.
*/
def makeGenomeBrowserFilesChannelSetupTss( inPaths, sampleSortedNames, sampleGenomeMap ) {
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
   }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inBedName = aSample + '-' + sampleGenomeMap[aSample]['name'] + '.tss_file.sorted.bed.gz'
        def tuple = new Tuple( fileMap[inBedName], [ 'sample': aSample, 'genome': sampleGenomeMap[aSample]['name'] ] )
        outTuples.add( tuple )
    }

    return( outTuples )
}


/*
** Set up a channel for getting merged peak browser file.
*/
def makeGenomeBrowserFilesChannelSetupMergedPeaks( inPaths, sampleSortedNames ) {
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inMergedPeaks = aSample + '-merged_peaks.bed'
        def tuple = new Tuple( fileMap[inMergedPeaks], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    return( outTuples )
}


/*
** Set up a channel for getting transposition site browser file.
*/
def makeGenomeBrowserFilesChannelSetupTranspositionSites( inPaths, sampleSortedNames ) {
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
    inPaths.each { aPath ->
        filesFound.add( aPath.getFileName().toString() )
    }
    filesExpected.each { aFile ->
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPath ->
        def aFile = aPath.getFileName().toString()
        def aSample = aFile.split( '-' )[0]
        if( aFile =~ /transposition_sites[.]bed[.]gz$/ ) {
            pathMap[aSample]['bed'] = aPath
        } else if( aFile =~ /transposition_sites[.]bed[.]gz[.]tbi$/ ) {
            pathMap[aSample]['tbi'] = aPath
        } else {
            println "Warning: makeGenomeBrowserFilesChannelSetupTranspositionSites: unexpected file \'${fileName}\'"
        }
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def tuple = new Tuple( pathMap[aSample]['bed'], pathMap[aSample]['tbi'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    return( outTuples )
}


/*
** Set up channel for making read alignment bigwig file.
*/
def makeGenomeBrowserFilesChannelSetupChromosomeSizes( inPaths, sampleSortedNames, sampleGenomeMap ) {
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def aGenomeJsonMap = sampleGenomeMap[aSample]
        def inGenomicIntervals = aSample + '-' + aGenomeJsonMap['name'] + '.chromosome_sizes.sorted.txt'
        def tuple = new Tuple( fileMap[inGenomicIntervals], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    return( outTuples )
}


/*
** Set up a channel for getting fragments.
*/
def getBandingScoresChannelSetupFragments( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected fragments.txt.gz files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-fragments.txt.gz'
        filesExpected.add( fileName )
    }
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-fragments.txt.gz.tbi'
        filesExpected.add( fileName )
    }
    def filesFound = []
    inPaths.each { aPath ->
        filesFound.add( aPath.getFileName().toString() )
    }
    filesExpected.each { aFile ->
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def pathMap = [:]
    sampleSortedNames.each { aSample ->
        pathMap[aSample] = [:]
    }
    inPaths.each { aPath ->
        def aFile = aPath.getFileName().toString()
        def aSample = aFile.split( '-' )[0]
        if( aFile =~ /fragments[.]txt[.]gz$/ ) {
            pathMap[aSample]['fragment'] = aPath
        } else if( aFile =~ /fragments[.]txt[.]gz[.]tbi$/ ) {
            pathMap[aSample]['fragmentTbi'] = aPath
        } else {
            println "Warning: getBandingScoresChannelSetupFragments: unexpected file \'${fileName}\'"
        }
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def tuple = new Tuple( pathMap[aSample]['fragment'], pathMap[aSample]['fragmentTbi'], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    return( outTuples )
}


/*
** Set up channel for making cell whitelist.
*/
def getBandingScoresChannelSetupCellWhitelist( inPaths, sampleSortedNames ) {
    /*
    ** Check for expected cell_whitelist.txt files.
    */
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        def fileName = aSample + '-called_cells_whitelist.txt'
        filesExpected.add( fileName )
    }
    def fileMap = getFileMap( inPaths )
    def filesFound = fileMap.keySet()
    filesExpected.each { aFile ->
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inCellWhitelist = aSample + '-called_cells_whitelist.txt'
        def tuple = new Tuple( fileMap[inCellWhitelist], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    return( outTuples )
}


def callMotifsChannelSetupMergedPeaks( inLists, sampleSortedNames, sampleGenomeMap ) {
    /*
    ** Make tuples of file names, sample names, and fasta and motifs file names.
    */
    def outTuples = []
	inLists.each { aList ->
		def aFile = aList[0].getFileName().toString()
		def aSample = aList[0].getFileName().toString().split( '-' )[0]
		def inMergedPeaks = aList[0]
		def gc_bin = aList[1]
		def aGenomeJsonMap = sampleGenomeMap[aSample]
        def fasta_exists = aGenomeJsonMap.containsKey( 'fasta' )
        def motifs_exists = aGenomeJsonMap.containsKey( 'motifs' )
        def fasta = fasta_exists ? aGenomeJsonMap['fasta'] : null
        if( fasta && !checkFile( fasta ) ) {
        	System.exit( -1 )
        }
        def motifs = motifs_exists ? aGenomeJsonMap['motifs'] : null
		if( motifs && !checkFile( motifs ) ) {
        	System.exit( -1 )
        }
        def tuple = new Tuple( inMergedPeaks, [ 'sample': aSample, 'genome': aGenomeJsonMap['name'], 'fasta': fasta, 'motifs': motifs, 'gc_bin': gc_bin ] )
        outTuples.add( tuple )	}
	
    return( outTuples )

}


def makeMotifMatrixChannelSetupPeakCalls( inPaths, sampleSortedNames, sampleGenomeMap ) {
	/*
	** Gather *.bb files into map keyed by sample name.
	*/
	def pathSetMap = [:]
	sampleSortedNames.each { aSample ->
		pathSetMap[aSample] = []
	}
	inPaths.each { aPath ->
		def aFile = aPath.getFileName().toString()
		def aSample = aPath.getFileName().toString().split( '-' )[0]
		pathSetMap[aSample].add( aPath )
	}

	/*
	** Check that each sample has the same number of files.
	** (A rudimentary test.)
	*/
	def countFiles = 0
	sampleSortedNames.each { aSample ->
		def fastaFlag = false
		def motifsFlag = false
		def aGenomeJsonMap = sampleGenomeMap[aSample]
        def fasta_exists = aGenomeJsonMap.containsKey( 'fasta' )
        def motifs_exists = aGenomeJsonMap.containsKey( 'motifs' )
        def fasta = fasta_exists ? aGenomeJsonMap['fasta'] : null
        if( fasta && checkFile( fasta ) ) {
        	fastaFlag = true
        }
        def motifs = motifs_exists ? aGenomeJsonMap['motifs'] : null
		if( motifs && checkFile( motifs ) ) {
        	motifsFlag = true
        }
        if( fastaFlag && motifsFlag ) {
			if( countFiles == 0 ) {
				countFiles = pathSetMap[aSample].size()
			}
			else {
				if( pathSetMap[aSample].size() != countFiles ) {
					printErr( "Error: inconsistent number of peaks calls files (by GC bin)" )
					System.exit( -1 )
				}
			}
		}
	}
	
	def outTuples = []
	sampleSortedNames.each { aSample ->
		pathSetMap[aSample] = pathSetMap[aSample].sort()
		def aGenomeJsonMap = sampleGenomeMap[aSample]
        def fasta_exists = aGenomeJsonMap.containsKey( 'fasta' )
        def motifs_exists = aGenomeJsonMap.containsKey( 'motifs' )
        def fasta = fasta_exists ? aGenomeJsonMap['fasta'] : null
        if( fasta && !checkFile( fasta ) ) {
        	System.exit( -1 )
        }
        def motifs = motifs_exists ? aGenomeJsonMap['motifs'] : null
		if( motifs && !checkFile( motifs ) ) {
        	System.exit( -1 )
        }
        def tuple = new Tuple( pathSetMap[aSample], [ 'sample': aSample, 'genome': aGenomeJsonMap['name'], 'fasta': fasta, 'motifs': motifs ] )
        outTuples.add( tuple )
	}

	return( outTuples )
}


def makeMotifMatrixChannelSetupMergedPeaks( inPaths, sampleSortedNames, sampleGenomeMap ) {
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
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inMergedPeaks = aSample + '-merged_peaks.bed'
		def aGenomeJsonMap = sampleGenomeMap[aSample]
        def fasta_exists = aGenomeJsonMap.containsKey( 'fasta' )
        def motifs_exists = aGenomeJsonMap.containsKey( 'motifs' )
        def fasta = fasta_exists ? aGenomeJsonMap['fasta'] : null
        if( fasta && !checkFile( fasta ) ) {
        	System.exit( -1 )
        }
        def motifs = motifs_exists ? aGenomeJsonMap['motifs'] : null
		if( motifs && !checkFile( motifs ) ) {
        	System.exit( -1 )
        }
        def tuple = new Tuple( fileMap[inMergedPeaks], [ 'sample': aSample, 'genome': aGenomeJsonMap['name'], 'fasta': fasta, 'motifs': motifs ] )
        outTuples.add( tuple )
    }

    return( outTuples )
}


def makeReducedDimensionMatrixChannelSetupPeakMatrix( inTuples, sampleSortedNames, sampleGenomeMap ) {
    /*
    ** Check for expected files.
    */
    def fileName
    def filesExpected = []
    sampleSortedNames.each { aSample ->
        fileName = aSample + '-peak_matrix.mtx.gz'
        filesExpected.add( fileName )
        fileName = aSample + '-peak_matrix.columns.txt'
        filesExpected.add( fileName )
        fileName = aSample + '-peak_matrix.rows.txt'
        filesExpected.add( fileName )
    }
    def filesFound = []
    inTuples.each { aTuple ->
    	aTuple.each { aPath ->
    		def aFile = aPath.getFileName().toString()
    		filesFound.add( aFile )
    	}
    }
    filesExpected.each { aFile ->
		if( !( aFile in filesFound ) ) {
			printErr( "Error: missing expected file \'${aFile}\' in channel" )
			System.exit( -1 )
		}
    }
    
	/*
	** Gather paths by sample.
	*/
    def pathSetMap = [:]
	sampleSortedNames.each { aSample ->
		pathSetMap[aSample] = []
    }
    inTuples.each { aTuple ->
    	aTuple.each { aPath ->
			def aFile = aPath.getFileName().toString()
			def aSample = aPath.getFileName().toString().split( '-' )[0]
			pathSetMap[aSample].add( aPath )    		
    	}
   	}
   	
	/*
	** Prepare output tuples.
	*/
	outTuples = []
    sampleSortedNames.each { aSample ->
		def aGenomeJsonMap = sampleGenomeMap[aSample]
        def blacklist_regions_exists = aGenomeJsonMap.containsKey( 'blacklist_regions' )
        def aBlackListRegion = blacklist_regions_exists ? aGenomeJsonMap['blacklist_regions'] : ''
    	tuple = new Tuple( pathSetMap[aSample], [ 'sample': aSample, 'genome': aGenomeJsonMap['name'], 'blacklist_regions': aBlackListRegion ] )
    	outTuples.add( tuple )
    }
	
	return( outTuples )
}


def makeReducedDimensionMatrixChannelSetupCountReport( inPaths, sampleSortedNames ) {
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
        if( !( aFile in filesFound ) ) {
            printErr( "Error: missing expected file \'${aFile}\' in channel" )
            System.exit( -1 )
        }
    }

    def outTuples = []
    sampleSortedNames.each { aSample ->
        def inCountReport = aSample + '-count_report.txt'
        def tuple = new Tuple( fileMap[inCountReport], [ 'sample': aSample ] )
        outTuples.add( tuple )
    }

    return( outTuples )
}


def experimentDashboardProcessChannelSetup( inPaths, sampleSortedNames, sampleGenomeMap ) {
	return( inPaths )
}


def makeMergedPlotFilesProcessChannelSetupCallCellsSummaryStats( inPaths, sampleSortedNames, sampleGenomeMap ) {
    return( inPaths )
}


def makeMergedPlotFilesProcessChannelSetupMakeMergedUmapPlots( inPaths, sampleSortedNames, sampleGenomeMap ) {
    return( inPaths )
}



