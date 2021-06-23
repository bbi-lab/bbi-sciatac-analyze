const log_data = {
"MM.Brain" :  `
***** PIPELINE READ STATS *****:

             Process       Starting reads         Ending reads               % lost      % of total lost
========================================================================================================
            Trimming             11633376             10602268                 8.86                 8.86
           Alignment             10602268              6922035                34.71                31.64
           Filtering              6922035              6006456                13.23                 7.87
       Deduplication              6006456              5328544                11.29                 5.83
     Gene assignment              5328544              4118836                22.70                10.40

Alignment details:
                                         Count              Percent
========================================================================================================
   Total reads processed:             10602268
   Reads uniquely mapped:              6006456                56.65
      Reads multi-mapped:               915579                 8.64
         Reads too short:              3330333                31.41
`,
"Sentinel" :  `
***** PIPELINE READ STATS *****:

             Process       Starting reads         Ending reads               % lost      % of total lost
========================================================================================================
            Trimming             14540997             13096032                 9.94                 9.94
           Alignment             13096032              8507325                35.04                31.56
           Filtering              8507325              6648851                21.85                12.78
       Deduplication              6648851              6049542                 9.01                 4.12
     Gene assignment              6049542              3946302                34.77                14.46

Alignment details:
                                         Count              Percent
========================================================================================================
   Total reads processed:             13096032
   Reads uniquely mapped:              6648851                50.77
      Reads multi-mapped:              1858474                14.19
         Reads too short:              4222495                32.24
`,
"NSM65" :  `
***** PIPELINE READ STATS *****:

             Process       Starting reads         Ending reads               % lost      % of total lost
========================================================================================================
            Trimming              5213162              4828837                 7.37                 7.37
           Alignment              4828837              2276870                52.85                48.95
           Filtering              2276870              1894102                16.81                 7.34
       Deduplication              1894102              1697691                10.37                 3.77
     Gene assignment              1697691               914113                46.16                15.03

Alignment details:
                                         Count              Percent
========================================================================================================
   Total reads processed:              4828837
   Reads uniquely mapped:              1894102                39.22
      Reads multi-mapped:               382768                 7.93
         Reads too short:              2378844                49.26
`,
"Barnyard" :  `
***** PIPELINE READ STATS *****:

             Process       Starting reads         Ending reads               % lost      % of total lost
========================================================================================================
            Trimming              9056751              8565855                 5.42                 5.42
           Alignment              8565855              7632398                10.90                10.31
           Filtering              7632398              6884976                 9.79                 8.25
       Deduplication              6884976              6091023                11.53                 8.77
     Gene assignment              6091023              4339858                28.75                19.34

Alignment details:
                                         Count              Percent
========================================================================================================
   Total reads processed:              8565855
   Reads uniquely mapped:              6884976                80.38
      Reads multi-mapped:               747422                 8.73
         Reads too short:               598672                 6.99
`
}
const full_log_data = {
"MM.Brain" :  `
BBI bbi-sci Pipeline Log

Nextflow version: 20.01.0
Pipeline version: 2.0.2
Git Repository, Version, Commit ID, Session ID: https://github.com/bbi-lab/bbi-sci.git, master, daf24eb76a7927c2a7e3acb9235add93485096f9, f09594ff-4462-420b-8770-96f919ad70f6

Command:
nextflow run bbi-sci -c experiment.config -resume

***** PARAMETERS *****:

    params.run_dir:               /net/bbi/vol1/seq/nextseq/200616_NS500773_0395_AHCGGYBGXF
    params.output_dir:            /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a
    params.sample_sheet:          /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/SampleSheet.csv
    params.p7_rows:               A B G H
    params.p5_cols:               8 9 7 6
    params.demux_out:             /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/demux_out
    params.level:                 3
    params.max_cores:             16
    params.samples:               false
    params.star_file:             /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/star_file.txt
    params.gene_file:             /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/gene_file.txt
    params.umi_cutoff:            100
    params.rt_barcode_file:       default
    params.hash_list:             false
    params.max_wells_per_sample:  20


Run started at: Thu Jun 18 05:31:12 PDT 2020

***** BEGIN PIPELINE *****:

** Start process 'check_sample_sheet' at: Thu Jun 18 05:31:12 PDT 2020

    Process versions:
        Python 3.6.4

    Process command:
        check_sample_sheet.py
            --sample_sheet /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/SampleSheet.csv
            --star_file /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/star_file.txt
            --level 3 --rt_barcode_file default
            --max_wells_per_samp 20

** End process 'check_sample_sheet' at: Thu Jun 18 05:31:15 PDT 2020

** Start process 'trim_fastqs' for MM.Brain-L002.fastq.gz at: Thu Jun 18 05:31:45 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore MM.Brain-L002.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: MM.Brain-L002.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA MM.Brain-L002.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 98.59 s (35 us/read; 1.71 M reads/minute).

=== Summary ===

Total reads processed:               2,805,192
Reads with adapters:                 1,998,538 (71.2%)
Reads written (passing filters):     2,805,192 (100.0%)

Total basepairs processed:   280,519,200 bp
Quality-trimmed:               8,383,450 bp (3.0%)
Total written (filtered):    191,893,040 bp (68.4%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 1998538 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 35.9%
  G: 24.1%
  T: 39.7%
  none/other: 0.3%


RUN STATISTICS FOR INPUT FILE: MM.Brain-L002.fastq.gz
=============================================
2805192 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	239295 (8.5%)

** End process 'trim_fastqs' at: Thu Jun 18 05:33:49 PDT 2020

** Start process 'trim_fastqs' for MM.Brain-L001.fastq.gz at: Thu Jun 18 05:31:40 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore MM.Brain-L001.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: MM.Brain-L001.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA MM.Brain-L001.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 99.81 s (35 us/read; 1.73 M reads/minute).

=== Summary ===

Total reads processed:               2,885,026
Reads with adapters:                 2,056,265 (71.3%)
Reads written (passing filters):     2,885,026 (100.0%)

Total basepairs processed:   288,502,600 bp
Quality-trimmed:               9,463,916 bp (3.3%)
Total written (filtered):    196,022,858 bp (67.9%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 2056265 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 36.1%
  G: 23.7%
  T: 39.9%
  none/other: 0.4%


RUN STATISTICS FOR INPUT FILE: MM.Brain-L001.fastq.gz
=============================================
2885026 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	256411 (8.9%)

** End process 'trim_fastqs' at: Thu Jun 18 05:33:44 PDT 2020

** Start process 'trim_fastqs' for MM.Brain-L004.fastq.gz at: Thu Jun 18 05:31:39 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore MM.Brain-L004.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: MM.Brain-L004.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA MM.Brain-L004.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 102.63 s (34 us/read; 1.75 M reads/minute).

=== Summary ===

Total reads processed:               2,999,798
Reads with adapters:                 2,144,232 (71.5%)
Reads written (passing filters):     2,999,798 (100.0%)

Total basepairs processed:   299,979,800 bp
Quality-trimmed:               9,636,464 bp (3.2%)
Total written (filtered):    203,069,346 bp (67.7%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 2144232 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 36.0%
  G: 23.7%
  T: 39.9%
  none/other: 0.4%


RUN STATISTICS FOR INPUT FILE: MM.Brain-L004.fastq.gz
=============================================
2999798 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	271325 (9.0%)

** End process 'trim_fastqs' at: Thu Jun 18 05:33:46 PDT 2020

** Start process 'trim_fastqs' for MM.Brain-L003.fastq.gz at: Thu Jun 18 05:31:30 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore MM.Brain-L003.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: MM.Brain-L003.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA MM.Brain-L003.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 108.00 s (37 us/read; 1.64 M reads/minute).

=== Summary ===

Total reads processed:               2,943,360
Reads with adapters:                 2,102,283 (71.4%)
Reads written (passing filters):     2,943,360 (100.0%)

Total basepairs processed:   294,336,000 bp
Quality-trimmed:               9,372,761 bp (3.2%)
Total written (filtered):    199,658,083 bp (67.8%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 2102283 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 36.0%
  G: 23.7%
  T: 39.9%
  none/other: 0.4%


RUN STATISTICS FOR INPUT FILE: MM.Brain-L003.fastq.gz
=============================================
2943360 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	264077 (9.0%)

** End process 'trim_fastqs' at: Thu Jun 18 05:33:41 PDT 2020

** Start process 'align_reads' for MM.Brain-L002_trimmed.fq.gz at: Thu Jun 18 05:34:10 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/mouse/mouse_star
            --readFilesIn MM.Brain-L002_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/MM.Brain-L002 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:10.09.00
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
        GTF download date:   2020.04.23:10.27.11


    Process output:
                                 Started job on |	Jun 18 05:34:10
                             Started mapping on |	Jun 18 05:38:00
                                    Finished on |	Jun 18 05:38:34
       Mapping speed, Million of reads per hour |	271.68

                          Number of input reads |	2565897
                      Average input read length |	72
                                    UNIQUE READS:
                   Uniquely mapped reads number |	1453376
                        Uniquely mapped reads % |	56.64%
                          Average mapped length |	69.42
                       Number of splices: Total |	18395
            Number of splices: Annotated (sjdb) |	13981
                       Number of splices: GT/AG |	17861
                       Number of splices: GC/AG |	459
                       Number of splices: AT/AC |	32
               Number of splices: Non-canonical |	43
                      Mismatch rate per base, % |	1.10%
                         Deletion rate per base |	0.02%
                        Deletion average length |	1.54
                        Insertion rate per base |	0.02%
                       Insertion average length |	1.11
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	219500
             % of reads mapped to multiple loci |	8.55%
        Number of reads mapped to too many loci |	26867
             % of reads mapped to too many loci |	1.05%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	31.54%
                     % of reads unmapped: other |	2.21%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:38:34 PDT 2020

** Start process 'align_reads' for MM.Brain-L001_trimmed.fq.gz at: Thu Jun 18 05:34:05 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/mouse/mouse_star
            --readFilesIn MM.Brain-L001_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/MM.Brain-L001 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:10.09.00
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
        GTF download date:   2020.04.23:10.27.11


    Process output:
                                 Started job on |	Jun 18 05:34:05
                             Started mapping on |	Jun 18 05:57:01
                                    Finished on |	Jun 18 05:57:24
       Mapping speed, Million of reads per hour |	411.44

                          Number of input reads |	2628615
                      Average input read length |	72
                                    UNIQUE READS:
                   Uniquely mapped reads number |	1490280
                        Uniquely mapped reads % |	56.69%
                          Average mapped length |	69.22
                       Number of splices: Total |	18860
            Number of splices: Annotated (sjdb) |	14274
                       Number of splices: GT/AG |	18289
                       Number of splices: GC/AG |	504
                       Number of splices: AT/AC |	23
               Number of splices: Non-canonical |	44
                      Mismatch rate per base, % |	1.00%
                         Deletion rate per base |	0.02%
                        Deletion average length |	1.55
                        Insertion rate per base |	0.02%
                       Insertion average length |	1.12
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	227340
             % of reads mapped to multiple loci |	8.65%
        Number of reads mapped to too many loci |	27866
             % of reads mapped to too many loci |	1.06%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	31.35%
                     % of reads unmapped: other |	2.25%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:57:24 PDT 2020

** Start process 'align_reads' for MM.Brain-L004_trimmed.fq.gz at: Thu Jun 18 05:34:07 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/mouse/mouse_star
            --readFilesIn MM.Brain-L004_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/MM.Brain-L004 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:10.09.00
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
        GTF download date:   2020.04.23:10.27.11


    Process output:
                                 Started job on |	Jun 18 05:34:07
                             Started mapping on |	Jun 18 05:58:18
                                    Finished on |	Jun 18 05:58:42
       Mapping speed, Million of reads per hour |	409.27

                          Number of input reads |	2728473
                      Average input read length |	71
                                    UNIQUE READS:
                   Uniquely mapped reads number |	1544275
                        Uniquely mapped reads % |	56.60%
                          Average mapped length |	69.10
                       Number of splices: Total |	19498
            Number of splices: Annotated (sjdb) |	14701
                       Number of splices: GT/AG |	18888
                       Number of splices: GC/AG |	522
                       Number of splices: AT/AC |	43
               Number of splices: Non-canonical |	45
                      Mismatch rate per base, % |	1.02%
                         Deletion rate per base |	0.02%
                        Deletion average length |	1.55
                        Insertion rate per base |	0.02%
                       Insertion average length |	1.12
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	236377
             % of reads mapped to multiple loci |	8.66%
        Number of reads mapped to too many loci |	29268
             % of reads mapped to too many loci |	1.07%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	31.43%
                     % of reads unmapped: other |	2.24%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:58:42 PDT 2020

** Start process 'align_reads' for MM.Brain-L003_trimmed.fq.gz at: Thu Jun 18 05:33:54 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/mouse/mouse_star
            --readFilesIn MM.Brain-L003_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/MM.Brain-L003 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:10.09.00
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
        GTF download date:   2020.04.23:10.27.11


    Process output:
                                 Started job on |	Jun 18 05:33:54
                             Started mapping on |	Jun 18 05:58:19
                                    Finished on |	Jun 18 05:58:58
       Mapping speed, Million of reads per hour |	247.32

                          Number of input reads |	2679283
                      Average input read length |	72
                                    UNIQUE READS:
                   Uniquely mapped reads number |	1518525
                        Uniquely mapped reads % |	56.68%
                          Average mapped length |	69.16
                       Number of splices: Total |	19320
            Number of splices: Annotated (sjdb) |	14626
                       Number of splices: GT/AG |	18735
                       Number of splices: GC/AG |	522
                       Number of splices: AT/AC |	26
               Number of splices: Non-canonical |	37
                      Mismatch rate per base, % |	1.00%
                         Deletion rate per base |	0.02%
                        Deletion average length |	1.55
                        Insertion rate per base |	0.02%
                       Insertion average length |	1.11
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	232362
             % of reads mapped to multiple loci |	8.67%
        Number of reads mapped to too many loci |	28965
             % of reads mapped to too many loci |	1.08%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	31.33%
                     % of reads unmapped: other |	2.24%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:58:58 PDT 2020

** Start process 'sort_and_filter' for MM.Brain-L002Aligned.out.bam at: Thu Jun 18 05:38:41 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'MM.Brain-L002Aligned.out.bam'
            | samtools sort -@ 10 - > 'MM.Brain-L002.bam'

    Process stats:
        sort_and_filter starting reads: 1672876
        sort_and_filter ending reads  : 1453376

** End process 'sort_and_filter' at: Thu Jun 18 05:39:21 PDT 2020

** Start process 'sort_and_filter' for MM.Brain-L001Aligned.out.bam at: Thu Jun 18 05:57:28 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'MM.Brain-L001Aligned.out.bam'
            | samtools sort -@ 10 - > 'MM.Brain-L001.bam'

    Process stats:
        sort_and_filter starting reads: 1717620
        sort_and_filter ending reads  : 1490280

** End process 'sort_and_filter' at: Thu Jun 18 05:58:02 PDT 2020

** Start process 'sort_and_filter' for MM.Brain-L004Aligned.out.bam at: Thu Jun 18 05:58:48 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'MM.Brain-L004Aligned.out.bam'
            | samtools sort -@ 10 - > 'MM.Brain-L004.bam'

    Process stats:
        sort_and_filter starting reads: 1780652
        sort_and_filter ending reads  : 1544275

** End process 'sort_and_filter' at: Thu Jun 18 05:59:22 PDT 2020

** Start process 'sort_and_filter' for MM.Brain-L003Aligned.out.bam at: Thu Jun 18 05:59:04 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'MM.Brain-L003Aligned.out.bam'
            | samtools sort -@ 10 - > 'MM.Brain-L003.bam'

    Process stats:
        sort_and_filter starting reads: 1750887
        sort_and_filter ending reads  : 1518525

** End process 'sort_and_filter' at: Thu Jun 18 05:59:37 PDT 2020

** Start process 'merge_bams' at: Thu Jun 18 06:21:00 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools merge MM.Brain.bam MM.Brain-L002.bam MM.Brain-L001.bam MM.Brain-L004.bam MM.Brain-L003.bam

** End process 'merge_bams' at: Thu Jun 18 06:21:20 PDT 2020

** Start processes 'remove duplicates, assign_genes, merge_assignment' at: Thu Jun 18 06:21:24 PDT 2020

    Process versions:
        bedtools v2.26.0
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.
        bamtools 2.2.3
        Python 3.6.4

    Process command:
        mkdir split_bams
        bamtools split -in MM.Brain.bam -reference -stub split_bams/split

        rmdup.py --bam in_bam --output_bam out.bam

        samtools view -c out.bam > split_bam_umi_count.txt

        bedtools bamtobed -i out.bam -split
                | sort -k1,1 -k2,2n -k3,3n -S 5G
                > "in_bam.bed"

        bedtools map
            -a in_bam.bed
            -b exon_index
            -nonamecheck -s -f 0.95 -c 7 -o distinct -delim "|"
        | bedtools map
            -a - -b gene_index
            -nonamecheck -s -f 0.95 -c 4 -o distinct -delim "|"
        | sort -k4,4 -k2,2n -k3,3n -S 5G
        | datamash
            -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8
        | assign-reads-to-genes.py gene_index
        | awk $3 == "exonic" || $3 == "intronic" {{
                split($1, arr, "|")
                printf "%s_%s_%s	%s	%s\n", arr[3], arr[4], arr[5], $2, $3
        }}
        | sort -k2,2 -k1,1 -S 5G > in_bam.txt

        cat logfile > merge_assignment.log
        cat split_bed > key.bed
        sort -m -k1,1 -k2,2 split_gene_assign > key_ga.txt

        datamash -g 1,2 count 2 < key_ga.txt
        | gzip > key.gz

    Process stats:
        remove_dups starting reads: 6006456
        remove_dups ending reads  : 5328544


        Read assignments:
            exonic                  2252018
            intronic                1866818

** End processes 'remove duplicates, assign_genes, merge_assignment' at: Thu Jun 18 06:58:41 PDT 2020

** Start process 'count_umis_by_sample' at: Thu Jun 18 06:58:45 PDT 2020

    Process versions:
        Python 3.6.4

    Process command:
        tabulate_per_cell_counts.py
            --gene_assignment_files MM.Brain_ga.txt
            --all_counts_file MM.Brain.UMIs.per.cell.barcode.txt
            --intron_counts_file MM.Brain.UMIs.per.cell.barcode.intronic.txt

    Process stats:
        Total cells                            : 312520
        Total cells > 100 reads                : 6980
        Total cells > 1000 reads               : 193
        Total reads in cells with > 100 reads  : 1941951

** End process 'count_umis_by_sample' at: Thu Jun 18 06:58:57 PDT 2020

** Start process 'make_matrix' at: Thu Jun 18 06:59:04 PDT 2020

    Process command:
        make_matrix.py <(zcat MM.Brain.gz)
            --gene_annotation "/net/bbi/vol1/data/genomes_stage/mouse/mouse_rna//latest.gene.annotations"
            --key "MM.Brain"
        cat /net/bbi/vol1/data/genomes_stage/mouse/mouse_rna//latest.gene.annotations > "MM.Brain.gene_annotations.txt"

** End process 'make_matrix' at: Thu Jun 18 06:59:32 PDT 2020

** Start process 'make_cds' at: Thu Jun 18 06:59:39 PDT 2020

    Process versions:
        R version 3.6.1 (2019-07-05) -- "Action of the Toes"
            monocle3 version [1] ‘0.2.1.2’

    Process command:
        make_cds.R
            "MM.Brain.umi_counts.mtx"
            "MM.Brain.cell_annotations.txt"
            "MM.Brain.gene_annotations.txt"
            "/net/bbi/vol1/data/genomes_stage/mouse/mouse_rna//latest.genes.bed"
            "MM.Brain"
            "100"

** End process 'make_cds' at: Thu Jun 18 07:00:12 PDT 2020

** Start process 'apply_garnett' at: Thu Jun 18 07:00:17 PDT 2020

ok

** End process 'apply_garnett' at: Thu Jun 18 07:01:37 PDT 2020

** Start process 'run_scrublet' at: Thu Jun 18 07:01:42 PDT 2020

    Process versions:
        Python 3.6.4
            scrublet  0.2.1

    Process command:
        run_scrublet.py --key MM.Brain --mat MM.Brain_for_scrub.mtx

** End process 'run_scrublet' at: Thu Jun 18 07:02:25 PDT 2020

** Start processes to generate qc metrics and dashboard at: Thu Jun 18 07:02:25 PDT 2020


** End processes generate qc metrics and dashboard at: Thu Jun 18 09:15:19 PDT 2020

***** END PIPELINE *****:

***** PIPELINE READ STATS *****:

             Process       Starting reads         Ending reads               % lost      % of total lost
========================================================================================================
            Trimming             11633376             10602268                 8.86                 8.86
           Alignment             10602268              6922035                34.71                31.64
           Filtering              6922035              6006456                13.23                 7.87
       Deduplication              6006456              5328544                11.29                 5.83
     Gene assignment              5328544              4118836                22.70                10.40

Alignment details:
                                         Count              Percent
========================================================================================================
   Total reads processed:             10602268
   Reads uniquely mapped:              6006456                56.65
      Reads multi-mapped:               915579                 8.64
         Reads too short:              3330333                31.41
`,
"Sentinel" :  `
BBI bbi-sci Pipeline Log

Nextflow version: 20.01.0
Pipeline version: 2.0.2
Git Repository, Version, Commit ID, Session ID: https://github.com/bbi-lab/bbi-sci.git, master, daf24eb76a7927c2a7e3acb9235add93485096f9, f09594ff-4462-420b-8770-96f919ad70f6

Command:
nextflow run bbi-sci -c experiment.config -resume

***** PARAMETERS *****:

    params.run_dir:               /net/bbi/vol1/seq/nextseq/200616_NS500773_0395_AHCGGYBGXF
    params.output_dir:            /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a
    params.sample_sheet:          /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/SampleSheet.csv
    params.p7_rows:               A B G H
    params.p5_cols:               8 9 7 6
    params.demux_out:             /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/demux_out
    params.level:                 3
    params.max_cores:             16
    params.samples:               false
    params.star_file:             /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/star_file.txt
    params.gene_file:             /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/gene_file.txt
    params.umi_cutoff:            100
    params.rt_barcode_file:       default
    params.hash_list:             false
    params.max_wells_per_sample:  20


Run started at: Thu Jun 18 05:31:12 PDT 2020

***** BEGIN PIPELINE *****:

** Start process 'check_sample_sheet' at: Thu Jun 18 05:31:12 PDT 2020

    Process versions:
        Python 3.6.4

    Process command:
        check_sample_sheet.py
            --sample_sheet /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/SampleSheet.csv
            --star_file /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/star_file.txt
            --level 3 --rt_barcode_file default
            --max_wells_per_samp 20

** End process 'check_sample_sheet' at: Thu Jun 18 05:31:15 PDT 2020

** Start process 'trim_fastqs' for Sentinel-L001.fastq.gz at: Thu Jun 18 05:31:26 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore Sentinel-L001.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: Sentinel-L001.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA Sentinel-L001.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 129.61 s (36 us/read; 1.67 M reads/minute).

=== Summary ===

Total reads processed:               3,605,170
Reads with adapters:                 2,555,232 (70.9%)
Reads written (passing filters):     3,605,170 (100.0%)

Total basepairs processed:   360,517,000 bp
Quality-trimmed:              12,493,791 bp (3.5%)
Total written (filtered):    242,878,478 bp (67.4%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 2555232 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 41.5%
  G: 22.0%
  T: 36.2%
  none/other: 0.3%


RUN STATISTICS FOR INPUT FILE: Sentinel-L001.fastq.gz
=============================================
3605170 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	358544 (9.9%)

** End process 'trim_fastqs' at: Thu Jun 18 05:34:03 PDT 2020

** Start process 'trim_fastqs' for Sentinel-L003.fastq.gz at: Thu Jun 18 05:31:23 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore Sentinel-L003.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: Sentinel-L003.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA Sentinel-L003.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 133.49 s (36 us/read; 1.65 M reads/minute).

=== Summary ===

Total reads processed:               3,681,893
Reads with adapters:                 2,614,194 (71.0%)
Reads written (passing filters):     3,681,893 (100.0%)

Total basepairs processed:   368,189,300 bp
Quality-trimmed:              12,395,451 bp (3.4%)
Total written (filtered):    247,651,276 bp (67.3%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 2614194 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 41.4%
  G: 22.0%
  T: 36.3%
  none/other: 0.3%


RUN STATISTICS FOR INPUT FILE: Sentinel-L003.fastq.gz
=============================================
3681893 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	369234 (10.0%)

** End process 'trim_fastqs' at: Thu Jun 18 05:34:08 PDT 2020

** Start process 'trim_fastqs' for Sentinel-L002.fastq.gz at: Thu Jun 18 05:31:32 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore Sentinel-L002.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: Sentinel-L002.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA Sentinel-L002.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 123.11 s (35 us/read; 1.71 M reads/minute).

=== Summary ===

Total reads processed:               3,501,449
Reads with adapters:                 2,476,511 (70.7%)
Reads written (passing filters):     3,501,449 (100.0%)

Total basepairs processed:   350,144,900 bp
Quality-trimmed:              11,118,569 bp (3.2%)
Total written (filtered):    237,784,965 bp (67.9%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 2476511 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 41.1%
  G: 22.6%
  T: 36.0%
  none/other: 0.3%


RUN STATISTICS FOR INPUT FILE: Sentinel-L002.fastq.gz
=============================================
3501449 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	335038 (9.6%)

** End process 'trim_fastqs' at: Thu Jun 18 05:34:04 PDT 2020

** Start process 'trim_fastqs' for Sentinel-L004.fastq.gz at: Thu Jun 18 05:31:38 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore Sentinel-L004.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: Sentinel-L004.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA Sentinel-L004.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 131.10 s (35 us/read; 1.72 M reads/minute).

=== Summary ===

Total reads processed:               3,752,485
Reads with adapters:                 2,667,369 (71.1%)
Reads written (passing filters):     3,752,485 (100.0%)

Total basepairs processed:   375,248,500 bp
Quality-trimmed:              12,839,804 bp (3.4%)
Total written (filtered):    251,718,697 bp (67.1%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 2667369 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 41.3%
  G: 22.0%
  T: 36.3%
  none/other: 0.4%


RUN STATISTICS FOR INPUT FILE: Sentinel-L004.fastq.gz
=============================================
3752485 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	382149 (10.2%)

** End process 'trim_fastqs' at: Thu Jun 18 05:34:17 PDT 2020

** Start process 'align_reads' for Sentinel-L001_trimmed.fq.gz at: Thu Jun 18 05:37:26 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/human/human_star
            --readFilesIn Sentinel-L001_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/Sentinel-L001 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:09.07.59
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
        GTF download date:   2020.04.23:09.58.35


    Process output:
                                 Started job on |	Jun 18 05:37:27
                             Started mapping on |	Jun 18 05:47:55
                                    Finished on |	Jun 18 05:48:40
       Mapping speed, Million of reads per hour |	259.73

                          Number of input reads |	3246626
                      Average input read length |	72
                                    UNIQUE READS:
                   Uniquely mapped reads number |	1651191
                        Uniquely mapped reads % |	50.86%
                          Average mapped length |	69.47
                       Number of splices: Total |	57865
            Number of splices: Annotated (sjdb) |	44275
                       Number of splices: GT/AG |	54354
                       Number of splices: GC/AG |	3338
                       Number of splices: AT/AC |	55
               Number of splices: Non-canonical |	118
                      Mismatch rate per base, % |	1.38%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.50
                        Insertion rate per base |	0.03%
                       Insertion average length |	1.29
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	459793
             % of reads mapped to multiple loci |	14.16%
        Number of reads mapped to too many loci |	41006
             % of reads mapped to too many loci |	1.26%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	32.18%
                     % of reads unmapped: other |	1.54%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:48:41 PDT 2020

** Start process 'align_reads' for Sentinel-L003_trimmed.fq.gz at: Thu Jun 18 05:37:27 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/human/human_star
            --readFilesIn Sentinel-L003_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/Sentinel-L003 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:09.07.59
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
        GTF download date:   2020.04.23:09.58.35


    Process output:
                                 Started job on |	Jun 18 05:37:28
                             Started mapping on |	Jun 18 05:52:09
                                    Finished on |	Jun 18 05:52:50
       Mapping speed, Million of reads per hour |	290.87

                          Number of input reads |	3312659
                      Average input read length |	72
                                    UNIQUE READS:
                   Uniquely mapped reads number |	1682514
                        Uniquely mapped reads % |	50.79%
                          Average mapped length |	69.44
                       Number of splices: Total |	58975
            Number of splices: Annotated (sjdb) |	45081
                       Number of splices: GT/AG |	55271
                       Number of splices: GC/AG |	3501
                       Number of splices: AT/AC |	79
               Number of splices: Non-canonical |	124
                      Mismatch rate per base, % |	1.38%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.48
                        Insertion rate per base |	0.03%
                       Insertion average length |	1.29
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	470698
             % of reads mapped to multiple loci |	14.21%
        Number of reads mapped to too many loci |	42472
             % of reads mapped to too many loci |	1.28%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	32.18%
                     % of reads unmapped: other |	1.54%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:52:50 PDT 2020

** Start process 'align_reads' for Sentinel-L002_trimmed.fq.gz at: Thu Jun 18 05:45:32 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/human/human_star
            --readFilesIn Sentinel-L002_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/Sentinel-L002 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:09.07.59
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
        GTF download date:   2020.04.23:09.58.35


    Process output:
                                 Started job on |	Jun 18 05:45:35
                             Started mapping on |	Jun 18 05:52:09
                                    Finished on |	Jun 18 05:52:50
       Mapping speed, Million of reads per hour |	278.03

                          Number of input reads |	3166411
                      Average input read length |	72
                                    UNIQUE READS:
                   Uniquely mapped reads number |	1603525
                        Uniquely mapped reads % |	50.64%
                          Average mapped length |	69.75
                       Number of splices: Total |	56517
            Number of splices: Annotated (sjdb) |	43109
                       Number of splices: GT/AG |	53024
                       Number of splices: GC/AG |	3324
                       Number of splices: AT/AC |	55
               Number of splices: Non-canonical |	114
                      Mismatch rate per base, % |	1.47%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.49
                        Insertion rate per base |	0.03%
                       Insertion average length |	1.27
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	447638
             % of reads mapped to multiple loci |	14.14%
        Number of reads mapped to too many loci |	40350
             % of reads mapped to too many loci |	1.27%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	32.46%
                     % of reads unmapped: other |	1.49%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:52:50 PDT 2020

** Start process 'align_reads' for Sentinel-L004_trimmed.fq.gz at: Thu Jun 18 05:47:58 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/human/human_star
            --readFilesIn Sentinel-L004_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/Sentinel-L004 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:09.07.59
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
        GTF download date:   2020.04.23:09.58.35


    Process output:
                                 Started job on |	Jun 18 05:47:58
                             Started mapping on |	Jun 18 06:05:16
                                    Finished on |	Jun 18 06:05:58
       Mapping speed, Million of reads per hour |	288.89

                          Number of input reads |	3370336
                      Average input read length |	72
                                    UNIQUE READS:
                   Uniquely mapped reads number |	1711621
                        Uniquely mapped reads % |	50.78%
                          Average mapped length |	69.37
                       Number of splices: Total |	59543
            Number of splices: Annotated (sjdb) |	45540
                       Number of splices: GT/AG |	55838
                       Number of splices: GC/AG |	3499
                       Number of splices: AT/AC |	69
               Number of splices: Non-canonical |	137
                      Mismatch rate per base, % |	1.40%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.50
                        Insertion rate per base |	0.03%
                       Insertion average length |	1.29
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	480345
             % of reads mapped to multiple loci |	14.25%
        Number of reads mapped to too many loci |	42925
             % of reads mapped to too many loci |	1.27%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	32.16%
                     % of reads unmapped: other |	1.53%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 06:05:58 PDT 2020

** Start process 'sort_and_filter' for Sentinel-L001Aligned.out.bam at: Thu Jun 18 05:48:53 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'Sentinel-L001Aligned.out.bam'
            | samtools sort -@ 10 - > 'Sentinel-L001.bam'

    Process stats:
        sort_and_filter starting reads: 2110984
        sort_and_filter ending reads  : 1651191

** End process 'sort_and_filter' at: Thu Jun 18 05:49:33 PDT 2020

** Start process 'sort_and_filter' for Sentinel-L003Aligned.out.bam at: Thu Jun 18 05:52:57 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'Sentinel-L003Aligned.out.bam'
            | samtools sort -@ 10 - > 'Sentinel-L003.bam'

    Process stats:
        sort_and_filter starting reads: 2153212
        sort_and_filter ending reads  : 1682514

** End process 'sort_and_filter' at: Thu Jun 18 05:53:36 PDT 2020

** Start process 'sort_and_filter' for Sentinel-L002Aligned.out.bam at: Thu Jun 18 05:52:59 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'Sentinel-L002Aligned.out.bam'
            | samtools sort -@ 10 - > 'Sentinel-L002.bam'

    Process stats:
        sort_and_filter starting reads: 2051163
        sort_and_filter ending reads  : 1603525

** End process 'sort_and_filter' at: Thu Jun 18 05:53:41 PDT 2020

** Start process 'sort_and_filter' for Sentinel-L004Aligned.out.bam at: Thu Jun 18 06:06:03 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'Sentinel-L004Aligned.out.bam'
            | samtools sort -@ 10 - > 'Sentinel-L004.bam'

    Process stats:
        sort_and_filter starting reads: 2191966
        sort_and_filter ending reads  : 1711621

** End process 'sort_and_filter' at: Thu Jun 18 06:06:44 PDT 2020

** Start process 'merge_bams' at: Thu Jun 18 06:21:04 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools merge Sentinel.bam Sentinel-L001.bam Sentinel-L003.bam Sentinel-L002.bam Sentinel-L004.bam

** End process 'merge_bams' at: Thu Jun 18 06:21:25 PDT 2020

** Start processes 'remove duplicates, assign_genes, merge_assignment' at: Thu Jun 18 06:21:32 PDT 2020

    Process versions:
        bedtools v2.26.0
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.
        bamtools 2.2.3
        Python 3.6.4

    Process command:
        mkdir split_bams
        bamtools split -in Sentinel.bam -reference -stub split_bams/split

        rmdup.py --bam in_bam --output_bam out.bam

        samtools view -c out.bam > split_bam_umi_count.txt

        bedtools bamtobed -i out.bam -split
                | sort -k1,1 -k2,2n -k3,3n -S 5G
                > "in_bam.bed"

        bedtools map
            -a in_bam.bed
            -b exon_index
            -nonamecheck -s -f 0.95 -c 7 -o distinct -delim "|"
        | bedtools map
            -a - -b gene_index
            -nonamecheck -s -f 0.95 -c 4 -o distinct -delim "|"
        | sort -k4,4 -k2,2n -k3,3n -S 5G
        | datamash
            -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8
        | assign-reads-to-genes.py gene_index
        | awk $3 == "exonic" || $3 == "intronic" {{
                split($1, arr, "|")
                printf "%s_%s_%s	%s	%s\n", arr[3], arr[4], arr[5], $2, $3
        }}
        | sort -k2,2 -k1,1 -S 5G > in_bam.txt

        cat logfile > merge_assignment.log
        cat split_bed > key.bed
        sort -m -k1,1 -k2,2 split_gene_assign > key_ga.txt

        datamash -g 1,2 count 2 < key_ga.txt
        | gzip > key.gz

    Process stats:
        remove_dups starting reads: 6648851
        remove_dups ending reads  : 6049542


        Read assignments:
            exonic                  2356928
            intronic                1589374

** End processes 'remove duplicates, assign_genes, merge_assignment' at: Thu Jun 18 06:58:48 PDT 2020

** Start process 'count_umis_by_sample' at: Thu Jun 18 06:58:53 PDT 2020

    Process versions:
        Python 3.6.4

    Process command:
        tabulate_per_cell_counts.py
            --gene_assignment_files Sentinel_ga.txt
            --all_counts_file Sentinel.UMIs.per.cell.barcode.txt
            --intron_counts_file Sentinel.UMIs.per.cell.barcode.intronic.txt

    Process stats:
        Total cells                            : 344793
        Total cells > 100 reads                : 5687
        Total cells > 1000 reads               : 234
        Total reads in cells with > 100 reads  : 1861290

** End process 'count_umis_by_sample' at: Thu Jun 18 06:59:02 PDT 2020

** Start process 'make_matrix' at: Thu Jun 18 06:59:08 PDT 2020

    Process command:
        make_matrix.py <(zcat Sentinel.gz)
            --gene_annotation "/net/bbi/vol1/data/genomes_stage/human/human_rna//latest.gene.annotations"
            --key "Sentinel"
        cat /net/bbi/vol1/data/genomes_stage/human/human_rna//latest.gene.annotations > "Sentinel.gene_annotations.txt"

** End process 'make_matrix' at: Thu Jun 18 06:59:37 PDT 2020

** Start process 'make_cds' at: Thu Jun 18 06:59:51 PDT 2020

    Process versions:
        R version 3.6.1 (2019-07-05) -- "Action of the Toes"
            monocle3 version [1] ‘0.2.1.2’

    Process command:
        make_cds.R
            "Sentinel.umi_counts.mtx"
            "Sentinel.cell_annotations.txt"
            "Sentinel.gene_annotations.txt"
            "/net/bbi/vol1/data/genomes_stage/human/human_rna//latest.genes.bed"
            "Sentinel"
            "100"

** End process 'make_cds' at: Thu Jun 18 07:00:18 PDT 2020

** Start process 'apply_garnett' at: Thu Jun 18 07:00:23 PDT 2020

No Garnett classifier provided for this sample

** End process 'apply_garnett' at: Thu Jun 18 07:00:39 PDT 2020

** Start process 'run_scrublet' at: Thu Jun 18 07:00:44 PDT 2020

    Process versions:
        Python 3.6.4
            scrublet  0.2.1

    Process command:
        run_scrublet.py --key Sentinel --mat Sentinel_for_scrub.mtx

** End process 'run_scrublet' at: Thu Jun 18 07:01:46 PDT 2020

** Start processes to generate qc metrics and dashboard at: Thu Jun 18 07:01:46 PDT 2020


** End processes generate qc metrics and dashboard at: Thu Jun 18 09:15:19 PDT 2020

***** END PIPELINE *****:

***** PIPELINE READ STATS *****:

             Process       Starting reads         Ending reads               % lost      % of total lost
========================================================================================================
            Trimming             14540997             13096032                 9.94                 9.94
           Alignment             13096032              8507325                35.04                31.56
           Filtering              8507325              6648851                21.85                12.78
       Deduplication              6648851              6049542                 9.01                 4.12
     Gene assignment              6049542              3946302                34.77                14.46

Alignment details:
                                         Count              Percent
========================================================================================================
   Total reads processed:             13096032
   Reads uniquely mapped:              6648851                50.77
      Reads multi-mapped:              1858474                14.19
         Reads too short:              4222495                32.24
`,
"3T3" :  `
BBI bbi-sci Pipeline Log

Nextflow version: 20.01.0
Pipeline version: 2.0.2
Git Repository, Version, Commit ID, Session ID: https://github.com/bbi-lab/bbi-sci.git, master, daf24eb76a7927c2a7e3acb9235add93485096f9, f09594ff-4462-420b-8770-96f919ad70f6

Command:
nextflow run bbi-sci -c experiment.config -resume

***** PARAMETERS *****:

    params.run_dir:               /net/bbi/vol1/seq/nextseq/200616_NS500773_0395_AHCGGYBGXF
    params.output_dir:            /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a
    params.sample_sheet:          /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/SampleSheet.csv
    params.p7_rows:               A B G H
    params.p5_cols:               8 9 7 6
    params.demux_out:             /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/demux_out
    params.level:                 3
    params.max_cores:             16
    params.samples:               false
    params.star_file:             /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/star_file.txt
    params.gene_file:             /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/gene_file.txt
    params.umi_cutoff:            100
    params.rt_barcode_file:       default
    params.hash_list:             false
    params.max_wells_per_sample:  20


Run started at: Thu Jun 18 05:31:12 PDT 2020

***** BEGIN PIPELINE *****:

** Start process 'check_sample_sheet' at: Thu Jun 18 05:31:12 PDT 2020

    Process versions:
        Python 3.6.4

    Process command:
        check_sample_sheet.py
            --sample_sheet /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/SampleSheet.csv
            --star_file /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/star_file.txt
            --level 3 --rt_barcode_file default
            --max_wells_per_samp 20

** End process 'check_sample_sheet' at: Thu Jun 18 05:31:15 PDT 2020

** Start process 'trim_fastqs' for 3T3.fq.part2-L003.fastq.gz at: Thu Jun 18 05:31:47 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore 3T3.fq.part2-L003.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: 3T3.fq.part2-L003.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA 3T3.fq.part2-L003.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 107.96 s (35 us/read; 1.73 M reads/minute).

=== Summary ===

Total reads processed:               3,105,203
Reads with adapters:                 1,791,218 (57.7%)
Reads written (passing filters):     3,105,203 (100.0%)

Total basepairs processed:   310,520,300 bp
Quality-trimmed:               7,037,558 bp (2.3%)
Total written (filtered):    246,656,461 bp (79.4%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 1791218 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 34.4%
  G: 20.0%
  T: 45.5%
  none/other: 0.2%


RUN STATISTICS FOR INPUT FILE: 3T3.fq.part2-L003.fastq.gz
=============================================
3105203 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	135340 (4.4%)

** End process 'trim_fastqs' at: Thu Jun 18 05:33:59 PDT 2020

** Start process 'trim_fastqs' for 3T3.fq.part1-L001.fastq.gz at: Thu Jun 18 05:31:36 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore 3T3.fq.part1-L001.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: 3T3.fq.part1-L001.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA 3T3.fq.part1-L001.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 326.18 s (36 us/read; 1.69 M reads/minute).

=== Summary ===

Total reads processed:               9,172,974
Reads with adapters:                 5,937,238 (64.7%)
Reads written (passing filters):     9,172,974 (100.0%)

Total basepairs processed:   917,297,400 bp
Quality-trimmed:              23,788,760 bp (2.6%)
Total written (filtered):    701,860,878 bp (76.5%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 5937238 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 35.2%
  G: 18.8%
  T: 45.8%
  none/other: 0.2%


RUN STATISTICS FOR INPUT FILE: 3T3.fq.part1-L001.fastq.gz
=============================================
9172974 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	465342 (5.1%)

** End process 'trim_fastqs' at: Thu Jun 18 05:37:56 PDT 2020

** Start process 'trim_fastqs' for 3T3.fq.part2-L001.fastq.gz at: Thu Jun 18 05:31:44 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore 3T3.fq.part2-L001.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: 3T3.fq.part2-L001.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA 3T3.fq.part2-L001.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 106.37 s (35 us/read; 1.72 M reads/minute).

=== Summary ===

Total reads processed:               3,052,202
Reads with adapters:                 1,756,870 (57.6%)
Reads written (passing filters):     3,052,202 (100.0%)

Total basepairs processed:   305,220,200 bp
Quality-trimmed:               7,114,452 bp (2.3%)
Total written (filtered):    242,597,344 bp (79.5%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 1756870 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 34.5%
  G: 19.9%
  T: 45.5%
  none/other: 0.1%


RUN STATISTICS FOR INPUT FILE: 3T3.fq.part2-L001.fastq.gz
=============================================
3052202 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	132140 (4.3%)

** End process 'trim_fastqs' at: Thu Jun 18 05:33:54 PDT 2020

** Start process 'trim_fastqs' for 3T3.fq.part1-L002.fastq.gz at: Thu Jun 18 05:31:26 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore 3T3.fq.part1-L002.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: 3T3.fq.part1-L002.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA 3T3.fq.part1-L002.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 339.79 s (38 us/read; 1.58 M reads/minute).

=== Summary ===

Total reads processed:               8,962,845
Reads with adapters:                 5,826,599 (65.0%)
Reads written (passing filters):     8,962,845 (100.0%)

Total basepairs processed:   896,284,500 bp
Quality-trimmed:              21,457,181 bp (2.4%)
Total written (filtered):    688,102,712 bp (76.8%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 5826599 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 34.8%
  G: 19.3%
  T: 45.7%
  none/other: 0.2%


RUN STATISTICS FOR INPUT FILE: 3T3.fq.part1-L002.fastq.gz
=============================================
8962845 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	444779 (5.0%)

** End process 'trim_fastqs' at: Thu Jun 18 05:38:00 PDT 2020

** Start process 'trim_fastqs' for 3T3.fq.part1-L004.fastq.gz at: Thu Jun 18 05:31:49 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore 3T3.fq.part1-L004.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: 3T3.fq.part1-L004.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA 3T3.fq.part1-L004.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 334.66 s (35 us/read; 1.70 M reads/minute).

=== Summary ===

Total reads processed:               9,502,342
Reads with adapters:                 6,162,602 (64.9%)
Reads written (passing filters):     9,502,342 (100.0%)

Total basepairs processed:   950,234,200 bp
Quality-trimmed:              24,551,964 bp (2.6%)
Total written (filtered):    725,042,440 bp (76.3%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 6162602 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 35.1%
  G: 18.8%
  T: 46.0%
  none/other: 0.2%


RUN STATISTICS FOR INPUT FILE: 3T3.fq.part1-L004.fastq.gz
=============================================
9502342 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	493281 (5.2%)

** End process 'trim_fastqs' at: Thu Jun 18 05:38:25 PDT 2020

** Start process 'trim_fastqs' for 3T3.fq.part2-L004.fastq.gz at: Thu Jun 18 05:31:31 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore 3T3.fq.part2-L004.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: 3T3.fq.part2-L004.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA 3T3.fq.part2-L004.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 109.67 s (35 us/read; 1.73 M reads/minute).

=== Summary ===

Total reads processed:               3,155,229
Reads with adapters:                 1,825,688 (57.9%)
Reads written (passing filters):     3,155,229 (100.0%)

Total basepairs processed:   315,522,900 bp
Quality-trimmed:               7,382,249 bp (2.3%)
Total written (filtered):    250,058,181 bp (79.3%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 1825688 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 34.4%
  G: 19.9%
  T: 45.5%
  none/other: 0.2%


RUN STATISTICS FOR INPUT FILE: 3T3.fq.part2-L004.fastq.gz
=============================================
3155229 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	139935 (4.4%)

** End process 'trim_fastqs' at: Thu Jun 18 05:33:49 PDT 2020

** Start process 'trim_fastqs' for 3T3.fq.part2-L002.fastq.gz at: Thu Jun 18 05:31:28 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore 3T3.fq.part2-L002.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: 3T3.fq.part2-L002.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA 3T3.fq.part2-L002.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 111.02 s (37 us/read; 1.61 M reads/minute).

=== Summary ===

Total reads processed:               2,983,908
Reads with adapters:                 1,734,054 (58.1%)
Reads written (passing filters):     2,983,908 (100.0%)

Total basepairs processed:   298,390,800 bp
Quality-trimmed:               6,523,843 bp (2.2%)
Total written (filtered):    237,781,880 bp (79.7%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 1734054 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 34.1%
  G: 20.3%
  T: 45.5%
  none/other: 0.2%


RUN STATISTICS FOR INPUT FILE: 3T3.fq.part2-L002.fastq.gz
=============================================
2983908 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	126247 (4.2%)

** End process 'trim_fastqs' at: Thu Jun 18 05:33:45 PDT 2020

** Start process 'trim_fastqs' for 3T3.fq.part1-L003.fastq.gz at: Thu Jun 18 05:31:35 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore 3T3.fq.part1-L003.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: 3T3.fq.part1-L003.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA 3T3.fq.part1-L003.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 329.46 s (35 us/read; 1.70 M reads/minute).

=== Summary ===

Total reads processed:               9,340,835
Reads with adapters:                 6,052,227 (64.8%)
Reads written (passing filters):     9,340,835 (100.0%)

Total basepairs processed:   934,083,500 bp
Quality-trimmed:              23,681,553 bp (2.5%)
Total written (filtered):    713,967,252 bp (76.4%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 6052227 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 35.1%
  G: 18.9%
  T: 45.9%
  none/other: 0.2%


RUN STATISTICS FOR INPUT FILE: 3T3.fq.part1-L003.fastq.gz
=============================================
9340835 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	480487 (5.1%)

** End process 'trim_fastqs' at: Thu Jun 18 05:38:00 PDT 2020

** Start process 'align_reads' for 3T3.fq.part2-L003_trimmed.fq.gz at: Thu Jun 18 05:37:24 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/mouse/mouse_star
            --readFilesIn 3T3.fq.part2-L003_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/3T3.fq.part2-L003 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:10.09.00
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
        GTF download date:   2020.04.23:10.27.11


    Process output:
                                 Started job on |	Jun 18 05:37:24
                             Started mapping on |	Jun 18 05:44:35
                                    Finished on |	Jun 18 05:45:08
       Mapping speed, Million of reads per hour |	323.99

                          Number of input reads |	2969863
                      Average input read length |	81
                                    UNIQUE READS:
                   Uniquely mapped reads number |	1983695
                        Uniquely mapped reads % |	66.79%
                          Average mapped length |	79.63
                       Number of splices: Total |	26427
            Number of splices: Annotated (sjdb) |	23364
                       Number of splices: GT/AG |	25756
                       Number of splices: GC/AG |	592
                       Number of splices: AT/AC |	36
               Number of splices: Non-canonical |	43
                      Mismatch rate per base, % |	0.77%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.74
                        Insertion rate per base |	0.02%
                       Insertion average length |	1.19
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	200422
             % of reads mapped to multiple loci |	6.75%
        Number of reads mapped to too many loci |	26099
             % of reads mapped to too many loci |	0.88%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	23.78%
                     % of reads unmapped: other |	1.79%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:45:08 PDT 2020

** Start process 'align_reads' for 3T3.fq.part1-L001_trimmed.fq.gz at: Thu Jun 18 05:38:37 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/mouse/mouse_star
            --readFilesIn 3T3.fq.part1-L001_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/3T3.fq.part1-L001 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:10.09.00
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
        GTF download date:   2020.04.23:10.27.11


    Process output:
                                 Started job on |	Jun 18 05:38:37
                             Started mapping on |	Jun 18 05:42:26
                                    Finished on |	Jun 18 05:43:41
       Mapping speed, Million of reads per hour |	417.97

                          Number of input reads |	8707632
                      Average input read length |	78
                                    UNIQUE READS:
                   Uniquely mapped reads number |	6309045
                        Uniquely mapped reads % |	72.45%
                          Average mapped length |	79.71
                       Number of splices: Total |	84772
            Number of splices: Annotated (sjdb) |	73634
                       Number of splices: GT/AG |	82539
                       Number of splices: GC/AG |	1968
                       Number of splices: AT/AC |	126
               Number of splices: Non-canonical |	139
                      Mismatch rate per base, % |	0.88%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.71
                        Insertion rate per base |	0.02%
                       Insertion average length |	1.18
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	683525
             % of reads mapped to multiple loci |	7.85%
        Number of reads mapped to too many loci |	82244
             % of reads mapped to too many loci |	0.94%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	16.84%
                     % of reads unmapped: other |	1.92%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:43:41 PDT 2020

** Start process 'align_reads' for 3T3.fq.part2-L001_trimmed.fq.gz at: Thu Jun 18 05:37:23 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/mouse/mouse_star
            --readFilesIn 3T3.fq.part2-L001_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/3T3.fq.part2-L001 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:10.09.00
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
        GTF download date:   2020.04.23:10.27.11


    Process output:
                                 Started job on |	Jun 18 05:37:23
                             Started mapping on |	Jun 18 05:47:15
                                    Finished on |	Jun 18 05:47:54
       Mapping speed, Million of reads per hour |	269.54

                          Number of input reads |	2920062
                      Average input read length |	81
                                    UNIQUE READS:
                   Uniquely mapped reads number |	1951870
                        Uniquely mapped reads % |	66.84%
                          Average mapped length |	79.64
                       Number of splices: Total |	26208
            Number of splices: Annotated (sjdb) |	23244
                       Number of splices: GT/AG |	25537
                       Number of splices: GC/AG |	596
                       Number of splices: AT/AC |	29
               Number of splices: Non-canonical |	46
                      Mismatch rate per base, % |	0.77%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.72
                        Insertion rate per base |	0.02%
                       Insertion average length |	1.19
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	196938
             % of reads mapped to multiple loci |	6.74%
        Number of reads mapped to too many loci |	25490
             % of reads mapped to too many loci |	0.87%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	23.75%
                     % of reads unmapped: other |	1.79%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:47:54 PDT 2020

** Start process 'align_reads' for 3T3.fq.part1-L002_trimmed.fq.gz at: Thu Jun 18 05:43:44 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/mouse/mouse_star
            --readFilesIn 3T3.fq.part1-L002_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/3T3.fq.part1-L002 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:10.09.00
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
        GTF download date:   2020.04.23:10.27.11


    Process output:
                                 Started job on |	Jun 18 05:43:44
                             Started mapping on |	Jun 18 05:47:35
                                    Finished on |	Jun 18 05:48:53
       Mapping speed, Million of reads per hour |	393.14

                          Number of input reads |	8518066
                      Average input read length |	79
                                    UNIQUE READS:
                   Uniquely mapped reads number |	6152327
                        Uniquely mapped reads % |	72.23%
                          Average mapped length |	79.79
                       Number of splices: Total |	82078
            Number of splices: Annotated (sjdb) |	70909
                       Number of splices: GT/AG |	79858
                       Number of splices: GC/AG |	1955
                       Number of splices: AT/AC |	113
               Number of splices: Non-canonical |	152
                      Mismatch rate per base, % |	1.00%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.70
                        Insertion rate per base |	0.02%
                       Insertion average length |	1.19
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	660360
             % of reads mapped to multiple loci |	7.75%
        Number of reads mapped to too many loci |	80260
             % of reads mapped to too many loci |	0.94%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	17.20%
                     % of reads unmapped: other |	1.88%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:48:53 PDT 2020

** Start process 'align_reads' for 3T3.fq.part1-L004_trimmed.fq.gz at: Thu Jun 18 05:48:56 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/mouse/mouse_star
            --readFilesIn 3T3.fq.part1-L004_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/3T3.fq.part1-L004 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:10.09.00
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
        GTF download date:   2020.04.23:10.27.11


    Process output:
                                 Started job on |	Jun 18 05:48:56
                             Started mapping on |	Jun 18 05:52:44
                                    Finished on |	Jun 18 05:54:02
       Mapping speed, Million of reads per hour |	415.80

                          Number of input reads |	9009061
                      Average input read length |	78
                                    UNIQUE READS:
                   Uniquely mapped reads number |	6509210
                        Uniquely mapped reads % |	72.25%
                          Average mapped length |	79.61
                       Number of splices: Total |	87775
            Number of splices: Annotated (sjdb) |	76024
                       Number of splices: GT/AG |	85416
                       Number of splices: GC/AG |	2111
                       Number of splices: AT/AC |	105
               Number of splices: Non-canonical |	143
                      Mismatch rate per base, % |	0.91%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.70
                        Insertion rate per base |	0.02%
                       Insertion average length |	1.19
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	710669
             % of reads mapped to multiple loci |	7.89%
        Number of reads mapped to too many loci |	86061
             % of reads mapped to too many loci |	0.96%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	16.99%
                     % of reads unmapped: other |	1.92%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:54:02 PDT 2020

** Start process 'align_reads' for 3T3.fq.part2-L004_trimmed.fq.gz at: Thu Jun 18 05:34:07 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/mouse/mouse_star
            --readFilesIn 3T3.fq.part2-L004_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/3T3.fq.part2-L004 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:10.09.00
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
        GTF download date:   2020.04.23:10.27.11


    Process output:
                                 Started job on |	Jun 18 05:34:07
                             Started mapping on |	Jun 18 05:58:18
                                    Finished on |	Jun 18 05:58:43
       Mapping speed, Million of reads per hour |	434.20

                          Number of input reads |	3015294
                      Average input read length |	81
                                    UNIQUE READS:
                   Uniquely mapped reads number |	2012765
                        Uniquely mapped reads % |	66.75%
                          Average mapped length |	79.52
                       Number of splices: Total |	26809
            Number of splices: Annotated (sjdb) |	23546
                       Number of splices: GT/AG |	26107
                       Number of splices: GC/AG |	637
                       Number of splices: AT/AC |	22
               Number of splices: Non-canonical |	43
                      Mismatch rate per base, % |	0.79%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.72
                        Insertion rate per base |	0.02%
                       Insertion average length |	1.20
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	204430
             % of reads mapped to multiple loci |	6.78%
        Number of reads mapped to too many loci |	26292
             % of reads mapped to too many loci |	0.87%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	23.80%
                     % of reads unmapped: other |	1.80%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:58:43 PDT 2020

** Start process 'align_reads' for 3T3.fq.part2-L002_trimmed.fq.gz at: Thu Jun 18 05:34:05 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/mouse/mouse_star
            --readFilesIn 3T3.fq.part2-L002_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/3T3.fq.part2-L002 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:10.09.00
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
        GTF download date:   2020.04.23:10.27.11


    Process output:
                                 Started job on |	Jun 18 05:34:05
                             Started mapping on |	Jun 18 06:04:09
                                    Finished on |	Jun 18 06:04:34
       Mapping speed, Million of reads per hour |	411.50

                          Number of input reads |	2857661
                      Average input read length |	81
                                    UNIQUE READS:
                   Uniquely mapped reads number |	1906562
                        Uniquely mapped reads % |	66.72%
                          Average mapped length |	79.70
                       Number of splices: Total |	25271
            Number of splices: Annotated (sjdb) |	22176
                       Number of splices: GT/AG |	24601
                       Number of splices: GC/AG |	587
                       Number of splices: AT/AC |	29
               Number of splices: Non-canonical |	54
                      Mismatch rate per base, % |	0.88%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.75
                        Insertion rate per base |	0.02%
                       Insertion average length |	1.19
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	191375
             % of reads mapped to multiple loci |	6.70%
        Number of reads mapped to too many loci |	24326
             % of reads mapped to too many loci |	0.85%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	23.98%
                     % of reads unmapped: other |	1.75%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 06:04:34 PDT 2020

** Start process 'align_reads' for 3T3.fq.part1-L003_trimmed.fq.gz at: Thu Jun 18 05:49:12 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/mouse/mouse_star
            --readFilesIn 3T3.fq.part1-L003_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/3T3.fq.part1-L003 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:10.09.00
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
        GTF download date:   2020.04.23:10.27.11


    Process output:
                                 Started job on |	Jun 18 05:49:19
                             Started mapping on |	Jun 18 06:06:48
                                    Finished on |	Jun 18 06:07:58
       Mapping speed, Million of reads per hour |	455.68

                          Number of input reads |	8860348
                      Average input read length |	78
                                    UNIQUE READS:
                   Uniquely mapped reads number |	6413962
                        Uniquely mapped reads % |	72.39%
                          Average mapped length |	79.69
                       Number of splices: Total |	86937
            Number of splices: Annotated (sjdb) |	75446
                       Number of splices: GT/AG |	84559
                       Number of splices: GC/AG |	2120
                       Number of splices: AT/AC |	104
               Number of splices: Non-canonical |	154
                      Mismatch rate per base, % |	0.88%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.70
                        Insertion rate per base |	0.03%
                       Insertion average length |	1.19
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	697421
             % of reads mapped to multiple loci |	7.87%
        Number of reads mapped to too many loci |	84032
             % of reads mapped to too many loci |	0.95%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	16.88%
                     % of reads unmapped: other |	1.91%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 06:07:58 PDT 2020

** Start process 'sort_and_filter' for 3T3.fq.part2-L003Aligned.out.bam at: Thu Jun 18 05:45:15 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 '3T3.fq.part2-L003Aligned.out.bam'
            | samtools sort -@ 10 - > '3T3.fq.part2-L003.bam'

    Process stats:
        sort_and_filter starting reads: 2184117
        sort_and_filter ending reads  : 1983695

** End process 'sort_and_filter' at: Thu Jun 18 05:46:21 PDT 2020

** Start process 'sort_and_filter' for 3T3.fq.part1-L001Aligned.out.bam at: Thu Jun 18 05:43:46 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 '3T3.fq.part1-L001Aligned.out.bam'
            | samtools sort -@ 10 - > '3T3.fq.part1-L001.bam'

    Process stats:
        sort_and_filter starting reads: 6992570
        sort_and_filter ending reads  : 6309045

** End process 'sort_and_filter' at: Thu Jun 18 05:46:43 PDT 2020

** Start process 'sort_and_filter' for 3T3.fq.part2-L001Aligned.out.bam at: Thu Jun 18 05:48:03 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 '3T3.fq.part2-L001Aligned.out.bam'
            | samtools sort -@ 10 - > '3T3.fq.part2-L001.bam'

    Process stats:
        sort_and_filter starting reads: 2148808
        sort_and_filter ending reads  : 1951870

** End process 'sort_and_filter' at: Thu Jun 18 05:48:55 PDT 2020

** Start process 'sort_and_filter' for 3T3.fq.part1-L002Aligned.out.bam at: Thu Jun 18 05:48:57 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 '3T3.fq.part1-L002Aligned.out.bam'
            | samtools sort -@ 10 - > '3T3.fq.part1-L002.bam'

    Process stats:
        sort_and_filter starting reads: 6812687
        sort_and_filter ending reads  : 6152327

** End process 'sort_and_filter' at: Thu Jun 18 05:51:30 PDT 2020

** Start process 'sort_and_filter' for 3T3.fq.part1-L004Aligned.out.bam at: Thu Jun 18 05:54:08 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 '3T3.fq.part1-L004Aligned.out.bam'
            | samtools sort -@ 10 - > '3T3.fq.part1-L004.bam'

    Process stats:
        sort_and_filter starting reads: 7219879
        sort_and_filter ending reads  : 6509210

** End process 'sort_and_filter' at: Thu Jun 18 05:56:47 PDT 2020

** Start process 'sort_and_filter' for 3T3.fq.part2-L004Aligned.out.bam at: Thu Jun 18 05:58:49 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 '3T3.fq.part2-L004Aligned.out.bam'
            | samtools sort -@ 10 - > '3T3.fq.part2-L004.bam'

    Process stats:
        sort_and_filter starting reads: 2217195
        sort_and_filter ending reads  : 2012765

** End process 'sort_and_filter' at: Thu Jun 18 05:59:34 PDT 2020

** Start process 'sort_and_filter' for 3T3.fq.part2-L002Aligned.out.bam at: Thu Jun 18 06:04:38 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 '3T3.fq.part2-L002Aligned.out.bam'
            | samtools sort -@ 10 - > '3T3.fq.part2-L002.bam'

    Process stats:
        sort_and_filter starting reads: 2097937
        sort_and_filter ending reads  : 1906562

** End process 'sort_and_filter' at: Thu Jun 18 06:05:25 PDT 2020

** Start process 'sort_and_filter' for 3T3.fq.part1-L003Aligned.out.bam at: Thu Jun 18 06:08:04 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 '3T3.fq.part1-L003Aligned.out.bam'
            | samtools sort -@ 10 - > '3T3.fq.part1-L003.bam'

    Process stats:
        sort_and_filter starting reads: 7111383
        sort_and_filter ending reads  : 6413962

** End process 'sort_and_filter' at: Thu Jun 18 06:10:37 PDT 2020

** Start process 'merge_bams' at: Thu Jun 18 06:21:05 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools merge 3T3.bam 3T3.fq.part2-L003.bam 3T3.fq.part1-L001.bam 3T3.fq.part2-L001.bam 3T3.fq.part1-L002.bam 3T3.fq.part1-L004.bam 3T3.fq.part2-L004.bam 3T3.fq.part2-L002.bam 3T3.fq.part1-L003.bam

** End process 'merge_bams' at: Thu Jun 18 06:22:48 PDT 2020

** Start processes 'remove duplicates, assign_genes, merge_assignment' at: Thu Jun 18 06:22:55 PDT 2020

    Process versions:
        bedtools v2.26.0
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.
        bamtools 2.2.3
        Python 3.6.4

    Process command:
        mkdir split_bams
        bamtools split -in 3T3.bam -reference -stub split_bams/split

        rmdup.py --bam in_bam --output_bam out.bam

        samtools view -c out.bam > split_bam_umi_count.txt

        bedtools bamtobed -i out.bam -split
                | sort -k1,1 -k2,2n -k3,3n -S 5G
                > "in_bam.bed"

        bedtools map
            -a in_bam.bed
            -b exon_index
            -nonamecheck -s -f 0.95 -c 7 -o distinct -delim "|"
        | bedtools map
            -a - -b gene_index
            -nonamecheck -s -f 0.95 -c 4 -o distinct -delim "|"
        | sort -k4,4 -k2,2n -k3,3n -S 5G
        | datamash
            -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8
        | assign-reads-to-genes.py gene_index
        | awk $3 == "exonic" || $3 == "intronic" {{
                split($1, arr, "|")
                printf "%s_%s_%s	%s	%s\n", arr[3], arr[4], arr[5], $2, $3
        }}
        | sort -k2,2 -k1,1 -S 5G > in_bam.txt

        cat logfile > merge_assignment.log
        cat split_bed > key.bed
        sort -m -k1,1 -k2,2 split_gene_assign > key_ga.txt

        datamash -g 1,2 count 2 < key_ga.txt
        | gzip > key.gz

    Process stats:
        remove_dups starting reads: 33239436
        remove_dups ending reads  : 29348794


        Read assignments:
            intronic                9713673
            exonic                 14424859

** End processes 'remove duplicates, assign_genes, merge_assignment' at: Thu Jun 18 07:00:15 PDT 2020

** Start process 'count_umis_by_sample' at: Thu Jun 18 07:00:18 PDT 2020

    Process versions:
        Python 3.6.4

    Process command:
        tabulate_per_cell_counts.py
            --gene_assignment_files 3T3_ga.txt
            --all_counts_file 3T3.UMIs.per.cell.barcode.txt
            --intron_counts_file 3T3.UMIs.per.cell.barcode.intronic.txt

    Process stats:
        Total cells                            : 433977
        Total cells > 100 reads                : 44043
        Total cells > 1000 reads               : 4407
        Total reads in cells with > 100 reads  : 21270006

** End process 'count_umis_by_sample' at: Thu Jun 18 07:01:10 PDT 2020

** Start process 'make_matrix' at: Thu Jun 18 07:01:12 PDT 2020

    Process command:
        make_matrix.py <(zcat 3T3.gz)
            --gene_annotation "/net/bbi/vol1/data/genomes_stage/mouse/mouse_rna//latest.gene.annotations"
            --key "3T3"
        cat /net/bbi/vol1/data/genomes_stage/mouse/mouse_rna//latest.gene.annotations > "3T3.gene_annotations.txt"

** End process 'make_matrix' at: Thu Jun 18 07:02:58 PDT 2020

** Start process 'make_cds' at: Thu Jun 18 07:03:04 PDT 2020

    Process versions:
        R version 3.6.1 (2019-07-05) -- "Action of the Toes"
            monocle3 version [1] ‘0.2.1.2’

    Process command:
        make_cds.R
            "3T3.umi_counts.mtx"
            "3T3.cell_annotations.txt"
            "3T3.gene_annotations.txt"
            "/net/bbi/vol1/data/genomes_stage/mouse/mouse_rna//latest.genes.bed"
            "3T3"
            "100"

** End process 'make_cds' at: Thu Jun 18 07:03:41 PDT 2020

** Start process 'apply_garnett' at: Thu Jun 18 07:03:49 PDT 2020

No Garnett classifier provided for this sample

** End process 'apply_garnett' at: Thu Jun 18 07:04:02 PDT 2020

** Start process 'run_scrublet' at: Thu Jun 18 07:04:07 PDT 2020

    Process versions:
        Python 3.6.4
            scrublet  0.2.1

    Process command:
        run_scrublet.py --key 3T3 --mat 3T3_for_scrub.mtx

** End process 'run_scrublet' at: Thu Jun 18 07:58:48 PDT 2020

** Start processes to generate qc metrics and dashboard at: Thu Jun 18 07:58:48 PDT 2020


** End processes generate qc metrics and dashboard at: Thu Jun 18 09:15:19 PDT 2020

***** END PIPELINE *****:

***** PIPELINE READ STATS *****:

             Process       Starting reads         Ending reads               % lost      % of total lost
========================================================================================================
            Trimming             49275538             46857987                 4.91                 4.91
           Alignment             46857987             36784576                21.50                20.44
           Filtering             36784576             33239436                 9.64                 7.19
       Deduplication             33239436             29348794                11.70                 7.90
     Gene assignment             29348794             24138532                17.75                10.57

Alignment details:
                                         Count              Percent
========================================================================================================
   Total reads processed:             46857987
   Reads uniquely mapped:             33239436                70.94
      Reads multi-mapped:              3545140                 7.57
         Reads too short:              8760394                18.70
`,
"NSM65" :  `
BBI bbi-sci Pipeline Log

Nextflow version: 20.01.0
Pipeline version: 2.0.2
Git Repository, Version, Commit ID, Session ID: https://github.com/bbi-lab/bbi-sci.git, master, daf24eb76a7927c2a7e3acb9235add93485096f9, f09594ff-4462-420b-8770-96f919ad70f6

Command:
nextflow run bbi-sci -c experiment.config -resume

***** PARAMETERS *****:

    params.run_dir:               /net/bbi/vol1/seq/nextseq/200616_NS500773_0395_AHCGGYBGXF
    params.output_dir:            /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a
    params.sample_sheet:          /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/SampleSheet.csv
    params.p7_rows:               A B G H
    params.p5_cols:               8 9 7 6
    params.demux_out:             /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/demux_out
    params.level:                 3
    params.max_cores:             16
    params.samples:               false
    params.star_file:             /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/star_file.txt
    params.gene_file:             /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/gene_file.txt
    params.umi_cutoff:            100
    params.rt_barcode_file:       default
    params.hash_list:             false
    params.max_wells_per_sample:  20


Run started at: Thu Jun 18 05:31:12 PDT 2020

***** BEGIN PIPELINE *****:

** Start process 'check_sample_sheet' at: Thu Jun 18 05:31:12 PDT 2020

    Process versions:
        Python 3.6.4

    Process command:
        check_sample_sheet.py
            --sample_sheet /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/SampleSheet.csv
            --star_file /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/star_file.txt
            --level 3 --rt_barcode_file default
            --max_wells_per_samp 20

** End process 'check_sample_sheet' at: Thu Jun 18 05:31:15 PDT 2020

** Start process 'trim_fastqs' for NSM65-L002.fastq.gz at: Thu Jun 18 05:31:21 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore NSM65-L002.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: NSM65-L002.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA NSM65-L002.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 45.71 s (36 us/read; 1.66 M reads/minute).

=== Summary ===

Total reads processed:               1,262,167
Reads with adapters:                   710,998 (56.3%)
Reads written (passing filters):     1,262,167 (100.0%)

Total basepairs processed:   126,216,700 bp
Quality-trimmed:               3,649,969 bp (2.9%)
Total written (filtered):     96,006,207 bp (76.1%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 710998 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 44.6%
  G: 22.1%
  T: 33.0%
  none/other: 0.3%


RUN STATISTICS FOR INPUT FILE: NSM65-L002.fastq.gz
=============================================
1262167 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	90877 (7.2%)

** End process 'trim_fastqs' at: Thu Jun 18 05:32:24 PDT 2020

** Start process 'trim_fastqs' for NSM65-L004.fastq.gz at: Thu Jun 18 05:31:22 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore NSM65-L004.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: NSM65-L004.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA NSM65-L004.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 49.21 s (37 us/read; 1.63 M reads/minute).

=== Summary ===

Total reads processed:               1,339,501
Reads with adapters:                   755,203 (56.4%)
Reads written (passing filters):     1,339,501 (100.0%)

Total basepairs processed:   133,950,100 bp
Quality-trimmed:               4,149,428 bp (3.1%)
Total written (filtered):    101,187,368 bp (75.5%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 755203 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 44.9%
  G: 21.8%
  T: 33.0%
  none/other: 0.3%


RUN STATISTICS FOR INPUT FILE: NSM65-L004.fastq.gz
=============================================
1339501 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	100958 (7.5%)

** End process 'trim_fastqs' at: Thu Jun 18 05:32:27 PDT 2020

** Start process 'trim_fastqs' for NSM65-L001.fastq.gz at: Thu Jun 18 05:31:29 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore NSM65-L001.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: NSM65-L001.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA NSM65-L001.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 47.31 s (37 us/read; 1.64 M reads/minute).

=== Summary ===

Total reads processed:               1,291,766
Reads with adapters:                   723,455 (56.0%)
Reads written (passing filters):     1,291,766 (100.0%)

Total basepairs processed:   129,176,600 bp
Quality-trimmed:               3,975,250 bp (3.1%)
Total written (filtered):     97,965,039 bp (75.8%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 723455 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 45.2%
  G: 21.6%
  T: 32.9%
  none/other: 0.3%


RUN STATISTICS FOR INPUT FILE: NSM65-L001.fastq.gz
=============================================
1291766 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	94716 (7.3%)

** End process 'trim_fastqs' at: Thu Jun 18 05:32:31 PDT 2020

** Start process 'trim_fastqs' for NSM65-L003.fastq.gz at: Thu Jun 18 05:31:39 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore NSM65-L003.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: NSM65-L003.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA NSM65-L003.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 45.25 s (34 us/read; 1.75 M reads/minute).

=== Summary ===

Total reads processed:               1,319,728
Reads with adapters:                   742,661 (56.3%)
Reads written (passing filters):     1,319,728 (100.0%)

Total basepairs processed:   131,972,800 bp
Quality-trimmed:               3,976,822 bp (3.0%)
Total written (filtered):     99,934,423 bp (75.7%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 742661 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 45.1%
  G: 21.8%
  T: 32.8%
  none/other: 0.3%


RUN STATISTICS FOR INPUT FILE: NSM65-L003.fastq.gz
=============================================
1319728 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	97774 (7.4%)

** End process 'trim_fastqs' at: Thu Jun 18 05:32:40 PDT 2020

** Start process 'align_reads' for NSM65-L002_trimmed.fq.gz at: Thu Jun 18 05:32:39 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/macaque/macaque_star
            --readFilesIn NSM65-L002_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/NSM65-L002 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp://ftp.ensembl.org/pub/release-99/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:13.04.37
        Non REF sequences removed.

        Genome GTF file URL: ftp://ftp.ensembl.org/pub/release-99/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.99.gtf.gz
        GTF download date:   2020.04.23:13.18.27


    Process output:
                                 Started job on |	Jun 18 05:32:39
                             Started mapping on |	Jun 18 05:36:54
                                    Finished on |	Jun 18 05:37:19
       Mapping speed, Million of reads per hour |	168.67

                          Number of input reads |	1171290
                      Average input read length |	79
                                    UNIQUE READS:
                   Uniquely mapped reads number |	457999
                        Uniquely mapped reads % |	39.10%
                          Average mapped length |	72.35
                       Number of splices: Total |	8586
            Number of splices: Annotated (sjdb) |	3046
                       Number of splices: GT/AG |	7597
                       Number of splices: GC/AG |	931
                       Number of splices: AT/AC |	11
               Number of splices: Non-canonical |	47
                      Mismatch rate per base, % |	2.57%
                         Deletion rate per base |	0.06%
                        Deletion average length |	2.00
                        Insertion rate per base |	0.06%
                       Insertion average length |	1.62
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	92999
             % of reads mapped to multiple loci |	7.94%
        Number of reads mapped to too many loci |	15471
             % of reads mapped to too many loci |	1.32%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	49.42%
                     % of reads unmapped: other |	2.21%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:37:19 PDT 2020

** Start process 'align_reads' for NSM65-L004_trimmed.fq.gz at: Thu Jun 18 05:32:40 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/macaque/macaque_star
            --readFilesIn NSM65-L004_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/NSM65-L004 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp://ftp.ensembl.org/pub/release-99/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:13.04.37
        Non REF sequences removed.

        Genome GTF file URL: ftp://ftp.ensembl.org/pub/release-99/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.99.gtf.gz
        GTF download date:   2020.04.23:13.18.27


    Process output:
                                 Started job on |	Jun 18 05:32:40
                             Started mapping on |	Jun 18 05:36:54
                                    Finished on |	Jun 18 05:37:20
       Mapping speed, Million of reads per hour |	171.49

                          Number of input reads |	1238543
                      Average input read length |	79
                                    UNIQUE READS:
                   Uniquely mapped reads number |	485594
                        Uniquely mapped reads % |	39.21%
                          Average mapped length |	72.13
                       Number of splices: Total |	9189
            Number of splices: Annotated (sjdb) |	3298
                       Number of splices: GT/AG |	8149
                       Number of splices: GC/AG |	990
                       Number of splices: AT/AC |	5
               Number of splices: Non-canonical |	45
                      Mismatch rate per base, % |	2.51%
                         Deletion rate per base |	0.06%
                        Deletion average length |	1.98
                        Insertion rate per base |	0.06%
                       Insertion average length |	1.62
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	98860
             % of reads mapped to multiple loci |	7.98%
        Number of reads mapped to too many loci |	16385
             % of reads mapped to too many loci |	1.32%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	49.19%
                     % of reads unmapped: other |	2.29%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:37:20 PDT 2020

** Start process 'align_reads' for NSM65-L001_trimmed.fq.gz at: Thu Jun 18 05:32:41 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/macaque/macaque_star
            --readFilesIn NSM65-L001_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/NSM65-L001 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp://ftp.ensembl.org/pub/release-99/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:13.04.37
        Non REF sequences removed.

        Genome GTF file URL: ftp://ftp.ensembl.org/pub/release-99/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.99.gtf.gz
        GTF download date:   2020.04.23:13.18.27


    Process output:
                                 Started job on |	Jun 18 05:32:41
                             Started mapping on |	Jun 18 06:01:38
                                    Finished on |	Jun 18 06:01:57
       Mapping speed, Million of reads per hour |	226.81

                          Number of input reads |	1197050
                      Average input read length |	79
                                    UNIQUE READS:
                   Uniquely mapped reads number |	470236
                        Uniquely mapped reads % |	39.28%
                          Average mapped length |	72.26
                       Number of splices: Total |	8889
            Number of splices: Annotated (sjdb) |	3236
                       Number of splices: GT/AG |	7847
                       Number of splices: GC/AG |	987
                       Number of splices: AT/AC |	6
               Number of splices: Non-canonical |	49
                      Mismatch rate per base, % |	2.49%
                         Deletion rate per base |	0.07%
                        Deletion average length |	2.02
                        Insertion rate per base |	0.06%
                       Insertion average length |	1.64
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	94435
             % of reads mapped to multiple loci |	7.89%
        Number of reads mapped to too many loci |	16126
             % of reads mapped to too many loci |	1.35%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	49.23%
                     % of reads unmapped: other |	2.25%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 06:01:57 PDT 2020

** Start process 'align_reads' for NSM65-L003_trimmed.fq.gz at: Thu Jun 18 05:32:51 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/macaque/macaque_star
            --readFilesIn NSM65-L003_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/NSM65-L003 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp://ftp.ensembl.org/pub/release-99/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:13.04.37
        Non REF sequences removed.

        Genome GTF file URL: ftp://ftp.ensembl.org/pub/release-99/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.99.gtf.gz
        GTF download date:   2020.04.23:13.18.27


    Process output:
                                 Started job on |	Jun 18 05:32:51
                             Started mapping on |	Jun 18 06:05:35
                                    Finished on |	Jun 18 06:05:54
       Mapping speed, Million of reads per hour |	231.53

                          Number of input reads |	1221954
                      Average input read length |	79
                                    UNIQUE READS:
                   Uniquely mapped reads number |	480273
                        Uniquely mapped reads % |	39.30%
                          Average mapped length |	72.18
                       Number of splices: Total |	9147
            Number of splices: Annotated (sjdb) |	3339
                       Number of splices: GT/AG |	8065
                       Number of splices: GC/AG |	1011
                       Number of splices: AT/AC |	13
               Number of splices: Non-canonical |	58
                      Mismatch rate per base, % |	2.48%
                         Deletion rate per base |	0.06%
                        Deletion average length |	2.00
                        Insertion rate per base |	0.06%
                       Insertion average length |	1.62
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	96474
             % of reads mapped to multiple loci |	7.90%
        Number of reads mapped to too many loci |	16260
             % of reads mapped to too many loci |	1.33%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	49.22%
                     % of reads unmapped: other |	2.25%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 06:05:54 PDT 2020

** Start process 'sort_and_filter' for NSM65-L002Aligned.out.bam at: Thu Jun 18 05:37:23 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'NSM65-L002Aligned.out.bam'
            | samtools sort -@ 10 - > 'NSM65-L002.bam'

    Process stats:
        sort_and_filter starting reads: 550998
        sort_and_filter ending reads  : 457999

** End process 'sort_and_filter' at: Thu Jun 18 05:37:36 PDT 2020

** Start process 'sort_and_filter' for NSM65-L004Aligned.out.bam at: Thu Jun 18 05:37:24 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'NSM65-L004Aligned.out.bam'
            | samtools sort -@ 10 - > 'NSM65-L004.bam'

    Process stats:
        sort_and_filter starting reads: 584454
        sort_and_filter ending reads  : 485594

** End process 'sort_and_filter' at: Thu Jun 18 05:37:39 PDT 2020

** Start process 'sort_and_filter' for NSM65-L001Aligned.out.bam at: Thu Jun 18 06:02:04 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'NSM65-L001Aligned.out.bam'
            | samtools sort -@ 10 - > 'NSM65-L001.bam'

    Process stats:
        sort_and_filter starting reads: 564671
        sort_and_filter ending reads  : 470236

** End process 'sort_and_filter' at: Thu Jun 18 06:02:18 PDT 2020

** Start process 'sort_and_filter' for NSM65-L003Aligned.out.bam at: Thu Jun 18 06:05:57 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'NSM65-L003Aligned.out.bam'
            | samtools sort -@ 10 - > 'NSM65-L003.bam'

    Process stats:
        sort_and_filter starting reads: 576747
        sort_and_filter ending reads  : 480273

** End process 'sort_and_filter' at: Thu Jun 18 06:06:09 PDT 2020

** Start process 'merge_bams' at: Thu Jun 18 06:20:59 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools merge NSM65.bam NSM65-L002.bam NSM65-L004.bam NSM65-L001.bam NSM65-L003.bam

** End process 'merge_bams' at: Thu Jun 18 06:21:06 PDT 2020

** Start processes 'remove duplicates, assign_genes, merge_assignment' at: Thu Jun 18 06:21:13 PDT 2020

    Process versions:
        bedtools v2.26.0
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.
        bamtools 2.2.3
        Python 3.6.4

    Process command:
        mkdir split_bams
        bamtools split -in NSM65.bam -reference -stub split_bams/split

        rmdup.py --bam in_bam --output_bam out.bam

        samtools view -c out.bam > split_bam_umi_count.txt

        bedtools bamtobed -i out.bam -split
                | sort -k1,1 -k2,2n -k3,3n -S 5G
                > "in_bam.bed"

        bedtools map
            -a in_bam.bed
            -b exon_index
            -nonamecheck -s -f 0.95 -c 7 -o distinct -delim "|"
        | bedtools map
            -a - -b gene_index
            -nonamecheck -s -f 0.95 -c 4 -o distinct -delim "|"
        | sort -k4,4 -k2,2n -k3,3n -S 5G
        | datamash
            -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8
        | assign-reads-to-genes.py gene_index
        | awk $3 == "exonic" || $3 == "intronic" {{
                split($1, arr, "|")
                printf "%s_%s_%s	%s	%s\n", arr[3], arr[4], arr[5], $2, $3
        }}
        | sort -k2,2 -k1,1 -S 5G > in_bam.txt

        cat logfile > merge_assignment.log
        cat split_bed > key.bed
        sort -m -k1,1 -k2,2 split_gene_assign > key_ga.txt

        datamash -g 1,2 count 2 < key_ga.txt
        | gzip > key.gz

    Process stats:
        remove_dups starting reads: 1894102
        remove_dups ending reads  : 1697691


        Read assignments:
            exonic                   148929
            intronic                 765184

** End processes 'remove duplicates, assign_genes, merge_assignment' at: Thu Jun 18 06:58:18 PDT 2020

** Start process 'count_umis_by_sample' at: Thu Jun 18 06:58:28 PDT 2020

    Process versions:
        Python 3.6.4

    Process command:
        tabulate_per_cell_counts.py
            --gene_assignment_files NSM65_ga.txt
            --all_counts_file NSM65.UMIs.per.cell.barcode.txt
            --intron_counts_file NSM65.UMIs.per.cell.barcode.intronic.txt

    Process stats:
        Total cells                            : 144833
        Total cells > 100 reads                : 1089
        Total cells > 1000 reads               : 32
        Total reads in cells with > 100 reads  : 289152

** End process 'count_umis_by_sample' at: Thu Jun 18 06:58:31 PDT 2020

** Start process 'make_matrix' at: Thu Jun 18 06:58:41 PDT 2020

    Process command:
        make_matrix.py <(zcat NSM65.gz)
            --gene_annotation "/net/bbi/vol1/data/genomes_stage/macaque/macaque_rna//latest.gene.annotations"
            --key "NSM65"
        cat /net/bbi/vol1/data/genomes_stage/macaque/macaque_rna//latest.gene.annotations > "NSM65.gene_annotations.txt"

** End process 'make_matrix' at: Thu Jun 18 06:58:49 PDT 2020

** Start process 'make_cds' at: Thu Jun 18 06:58:54 PDT 2020

    Process versions:
        R version 3.6.1 (2019-07-05) -- "Action of the Toes"
            monocle3 version [1] ‘0.2.1.2’

    Process command:
        make_cds.R
            "NSM65.umi_counts.mtx"
            "NSM65.cell_annotations.txt"
            "NSM65.gene_annotations.txt"
            "/net/bbi/vol1/data/genomes_stage/macaque/macaque_rna//latest.genes.bed"
            "NSM65"
            "100"

** End process 'make_cds' at: Thu Jun 18 06:59:34 PDT 2020

** Start process 'apply_garnett' at: Thu Jun 18 06:59:38 PDT 2020

No Garnett classifier provided for this sample

** End process 'apply_garnett' at: Thu Jun 18 06:59:55 PDT 2020

** Start process 'run_scrublet' at: Thu Jun 18 07:00:17 PDT 2020

    Process versions:
        Python 3.6.4
            scrublet  0.2.1

    Process command:
        run_scrublet.py --key NSM65 --mat NSM65_for_scrub.mtx

** End process 'run_scrublet' at: Thu Jun 18 07:01:29 PDT 2020

** Start processes to generate qc metrics and dashboard at: Thu Jun 18 07:01:29 PDT 2020


** End processes generate qc metrics and dashboard at: Thu Jun 18 09:15:19 PDT 2020

***** END PIPELINE *****:

***** PIPELINE READ STATS *****:

             Process       Starting reads         Ending reads               % lost      % of total lost
========================================================================================================
            Trimming              5213162              4828837                 7.37                 7.37
           Alignment              4828837              2276870                52.85                48.95
           Filtering              2276870              1894102                16.81                 7.34
       Deduplication              1894102              1697691                10.37                 3.77
     Gene assignment              1697691               914113                46.16                15.03

Alignment details:
                                         Count              Percent
========================================================================================================
   Total reads processed:              4828837
   Reads uniquely mapped:              1894102                39.22
      Reads multi-mapped:               382768                 7.93
         Reads too short:              2378844                49.26
`,
"Barnyard" :  `
BBI bbi-sci Pipeline Log

Nextflow version: 20.01.0
Pipeline version: 2.0.2
Git Repository, Version, Commit ID, Session ID: https://github.com/bbi-lab/bbi-sci.git, master, daf24eb76a7927c2a7e3acb9235add93485096f9, f09594ff-4462-420b-8770-96f919ad70f6

Command:
nextflow run bbi-sci -c experiment.config -resume

***** PARAMETERS *****:

    params.run_dir:               /net/bbi/vol1/seq/nextseq/200616_NS500773_0395_AHCGGYBGXF
    params.output_dir:            /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a
    params.sample_sheet:          /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/SampleSheet.csv
    params.p7_rows:               A B G H
    params.p5_cols:               8 9 7 6
    params.demux_out:             /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/demux_out
    params.level:                 3
    params.max_cores:             16
    params.samples:               false
    params.star_file:             /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/star_file.txt
    params.gene_file:             /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/gene_file.txt
    params.umi_cutoff:            100
    params.rt_barcode_file:       default
    params.hash_list:             false
    params.max_wells_per_sample:  20


Run started at: Thu Jun 18 05:31:12 PDT 2020

***** BEGIN PIPELINE *****:

** Start process 'check_sample_sheet' at: Thu Jun 18 05:31:12 PDT 2020

    Process versions:
        Python 3.6.4

    Process command:
        check_sample_sheet.py
            --sample_sheet /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/SampleSheet.csv
            --star_file /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/star_file.txt
            --level 3 --rt_barcode_file default
            --max_wells_per_samp 20

** End process 'check_sample_sheet' at: Thu Jun 18 05:31:15 PDT 2020

** Start process 'trim_fastqs' for Barnyard-L003.fastq.gz at: Thu Jun 18 05:31:34 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore Barnyard-L003.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: Barnyard-L003.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA Barnyard-L003.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 78.81 s (34 us/read; 1.74 M reads/minute).

=== Summary ===

Total reads processed:               2,286,192
Reads with adapters:                 1,491,408 (65.2%)
Reads written (passing filters):     2,286,192 (100.0%)

Total basepairs processed:   228,619,200 bp
Quality-trimmed:               6,143,636 bp (2.7%)
Total written (filtered):    173,215,679 bp (75.8%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 1491408 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 44.1%
  G: 17.0%
  T: 38.7%
  none/other: 0.1%


RUN STATISTICS FOR INPUT FILE: Barnyard-L003.fastq.gz
=============================================
2286192 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	124222 (5.4%)

** End process 'trim_fastqs' at: Thu Jun 18 05:33:14 PDT 2020

** Start process 'trim_fastqs' for Barnyard-L001.fastq.gz at: Thu Jun 18 05:31:24 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore Barnyard-L001.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: Barnyard-L001.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA Barnyard-L001.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 83.37 s (37 us/read; 1.62 M reads/minute).

=== Summary ===

Total reads processed:               2,244,860
Reads with adapters:                 1,464,113 (65.2%)
Reads written (passing filters):     2,244,860 (100.0%)

Total basepairs processed:   224,486,000 bp
Quality-trimmed:               6,205,433 bp (2.8%)
Total written (filtered):    170,119,996 bp (75.8%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 1464113 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 44.3%
  G: 16.9%
  T: 38.7%
  none/other: 0.1%


RUN STATISTICS FOR INPUT FILE: Barnyard-L001.fastq.gz
=============================================
2244860 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	121348 (5.4%)

** End process 'trim_fastqs' at: Thu Jun 18 05:33:09 PDT 2020

** Start process 'trim_fastqs' for Barnyard-L004.fastq.gz at: Thu Jun 18 05:31:48 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore Barnyard-L004.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: Barnyard-L004.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA Barnyard-L004.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 82.07 s (35 us/read; 1.70 M reads/minute).

=== Summary ===

Total reads processed:               2,328,716
Reads with adapters:                 1,521,951 (65.4%)
Reads written (passing filters):     2,328,716 (100.0%)

Total basepairs processed:   232,871,600 bp
Quality-trimmed:               6,369,646 bp (2.7%)
Total written (filtered):    176,022,365 bp (75.6%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 1521951 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 44.1%
  G: 17.0%
  T: 38.8%
  none/other: 0.1%


RUN STATISTICS FOR INPUT FILE: Barnyard-L004.fastq.gz
=============================================
2328716 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	128246 (5.5%)

** End process 'trim_fastqs' at: Thu Jun 18 05:33:30 PDT 2020

** Start process 'trim_fastqs' for Barnyard-L002.fastq.gz at: Thu Jun 18 05:31:50 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore Barnyard-L002.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: Barnyard-L002.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA Barnyard-L002.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 79.23 s (36 us/read; 1.66 M reads/minute).

=== Summary ===

Total reads processed:               2,196,983
Reads with adapters:                 1,445,610 (65.8%)
Reads written (passing filters):     2,196,983 (100.0%)

Total basepairs processed:   219,698,300 bp
Quality-trimmed:               5,557,120 bp (2.5%)
Total written (filtered):    166,905,169 bp (76.0%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 1445610 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 43.7%
  G: 17.4%
  T: 38.8%
  none/other: 0.2%


RUN STATISTICS FOR INPUT FILE: Barnyard-L002.fastq.gz
=============================================
2196983 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	117080 (5.3%)

** End process 'trim_fastqs' at: Thu Jun 18 05:33:29 PDT 2020

** Start process 'align_reads' for Barnyard-L003_trimmed.fq.gz at: Thu Jun 18 05:33:26 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/barnyard/barnyard_star
            --readFilesIn Barnyard-L003_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/Barnyard-L003 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:

        FASTA download date:
        Non REF sequences removed.


        GTF download date:


    Process output:
                                 Started job on |	Jun 18 05:33:26
                             Started mapping on |	Jun 18 06:16:19
                                    Finished on |	Jun 18 06:16:42
       Mapping speed, Million of reads per hour |	338.40

                          Number of input reads |	2161970
                      Average input read length |	78
                                    UNIQUE READS:
                   Uniquely mapped reads number |	1739712
                        Uniquely mapped reads % |	80.47%
                          Average mapped length |	80.20
                       Number of splices: Total |	20429
            Number of splices: Annotated (sjdb) |	14858
                       Number of splices: GT/AG |	18912
                       Number of splices: GC/AG |	1450
                       Number of splices: AT/AC |	23
               Number of splices: Non-canonical |	44
                      Mismatch rate per base, % |	0.90%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.55
                        Insertion rate per base |	0.03%
                       Insertion average length |	1.13
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	188212
             % of reads mapped to multiple loci |	8.71%
        Number of reads mapped to too many loci |	29852
             % of reads mapped to too many loci |	1.38%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	6.91%
                     % of reads unmapped: other |	2.54%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 06:16:42 PDT 2020

** Start process 'align_reads' for Barnyard-L001_trimmed.fq.gz at: Thu Jun 18 05:33:35 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/barnyard/barnyard_star
            --readFilesIn Barnyard-L001_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/Barnyard-L001 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:

        FASTA download date:
        Non REF sequences removed.


        GTF download date:


    Process output:
                                 Started job on |	Jun 18 05:33:37
                             Started mapping on |	Jun 18 06:18:15
                                    Finished on |	Jun 18 06:18:37
       Mapping speed, Million of reads per hour |	347.48

                          Number of input reads |	2123512
                      Average input read length |	78
                                    UNIQUE READS:
                   Uniquely mapped reads number |	1708516
                        Uniquely mapped reads % |	80.46%
                          Average mapped length |	80.19
                       Number of splices: Total |	20037
            Number of splices: Annotated (sjdb) |	14612
                       Number of splices: GT/AG |	18609
                       Number of splices: GC/AG |	1358
                       Number of splices: AT/AC |	35
               Number of splices: Non-canonical |	35
                      Mismatch rate per base, % |	0.90%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.54
                        Insertion rate per base |	0.03%
                       Insertion average length |	1.12
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	185243
             % of reads mapped to multiple loci |	8.72%
        Number of reads mapped to too many loci |	29499
             % of reads mapped to too many loci |	1.39%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	6.90%
                     % of reads unmapped: other |	2.53%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 06:18:37 PDT 2020

** Start process 'align_reads' for Barnyard-L004_trimmed.fq.gz at: Thu Jun 18 05:33:41 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/barnyard/barnyard_star
            --readFilesIn Barnyard-L004_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/Barnyard-L004 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:

        FASTA download date:
        Non REF sequences removed.


        GTF download date:


    Process output:
                                 Started job on |	Jun 18 05:33:41
                             Started mapping on |	Jun 18 06:18:15
                                    Finished on |	Jun 18 06:18:37
       Mapping speed, Million of reads per hour |	360.08

                          Number of input reads |	2200470
                      Average input read length |	78
                                    UNIQUE READS:
                   Uniquely mapped reads number |	1766277
                        Uniquely mapped reads % |	80.27%
                          Average mapped length |	80.07
                       Number of splices: Total |	20881
            Number of splices: Annotated (sjdb) |	15268
                       Number of splices: GT/AG |	19351
                       Number of splices: GC/AG |	1465
                       Number of splices: AT/AC |	28
               Number of splices: Non-canonical |	37
                      Mismatch rate per base, % |	0.93%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.54
                        Insertion rate per base |	0.03%
                       Insertion average length |	1.13
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	192922
             % of reads mapped to multiple loci |	8.77%
        Number of reads mapped to too many loci |	30589
             % of reads mapped to too many loci |	1.39%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	7.01%
                     % of reads unmapped: other |	2.56%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 06:18:37 PDT 2020

** Start process 'align_reads' for Barnyard-L002_trimmed.fq.gz at: Thu Jun 18 05:34:02 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/barnyard/barnyard_star
            --readFilesIn Barnyard-L002_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/Barnyard-L002 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:

        FASTA download date:
        Non REF sequences removed.


        GTF download date:


    Process output:
                                 Started job on |	Jun 18 05:34:03
                             Started mapping on |	Jun 18 06:19:44
                                    Finished on |	Jun 18 06:20:04
       Mapping speed, Million of reads per hour |	374.38

                          Number of input reads |	2079903
                      Average input read length |	78
                                    UNIQUE READS:
                   Uniquely mapped reads number |	1670471
                        Uniquely mapped reads % |	80.31%
                          Average mapped length |	80.19
                       Number of splices: Total |	19789
            Number of splices: Annotated (sjdb) |	14338
                       Number of splices: GT/AG |	18290
                       Number of splices: GC/AG |	1425
                       Number of splices: AT/AC |	28
               Number of splices: Non-canonical |	46
                      Mismatch rate per base, % |	1.02%
                         Deletion rate per base |	0.03%
                        Deletion average length |	1.55
                        Insertion rate per base |	0.03%
                       Insertion average length |	1.12
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	181045
             % of reads mapped to multiple loci |	8.70%
        Number of reads mapped to too many loci |	28318
             % of reads mapped to too many loci |	1.36%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	7.14%
                     % of reads unmapped: other |	2.48%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 06:20:04 PDT 2020

** Start process 'sort_and_filter' for Barnyard-L003Aligned.out.bam at: Thu Jun 18 06:16:48 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'Barnyard-L003Aligned.out.bam'
            | samtools sort -@ 10 - > 'Barnyard-L003.bam'

    Process stats:
        sort_and_filter starting reads: 1927924
        sort_and_filter ending reads  : 1739712

** End process 'sort_and_filter' at: Thu Jun 18 06:17:29 PDT 2020

** Start process 'sort_and_filter' for Barnyard-L001Aligned.out.bam at: Thu Jun 18 06:18:43 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'Barnyard-L001Aligned.out.bam'
            | samtools sort -@ 10 - > 'Barnyard-L001.bam'

    Process stats:
        sort_and_filter starting reads: 1893759
        sort_and_filter ending reads  : 1708516

** End process 'sort_and_filter' at: Thu Jun 18 06:19:24 PDT 2020

** Start process 'sort_and_filter' for Barnyard-L004Aligned.out.bam at: Thu Jun 18 06:18:44 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'Barnyard-L004Aligned.out.bam'
            | samtools sort -@ 10 - > 'Barnyard-L004.bam'

    Process stats:
        sort_and_filter starting reads: 1959199
        sort_and_filter ending reads  : 1766277

** End process 'sort_and_filter' at: Thu Jun 18 06:19:27 PDT 2020

** Start process 'sort_and_filter' for Barnyard-L002Aligned.out.bam at: Thu Jun 18 06:20:08 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 'Barnyard-L002Aligned.out.bam'
            | samtools sort -@ 10 - > 'Barnyard-L002.bam'

    Process stats:
        sort_and_filter starting reads: 1851516
        sort_and_filter ending reads  : 1670471

** End process 'sort_and_filter' at: Thu Jun 18 06:20:49 PDT 2020

** Start process 'merge_bams' at: Thu Jun 18 06:21:01 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools merge Barnyard.bam Barnyard-L003.bam Barnyard-L001.bam Barnyard-L004.bam Barnyard-L002.bam

** End process 'merge_bams' at: Thu Jun 18 06:21:25 PDT 2020

** Start processes 'remove duplicates, assign_genes, merge_assignment' at: Thu Jun 18 06:21:29 PDT 2020

    Process versions:
        bedtools v2.26.0
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.
        bamtools 2.2.3
        Python 3.6.4

    Process command:
        mkdir split_bams
        bamtools split -in Barnyard.bam -reference -stub split_bams/split

        rmdup.py --bam in_bam --output_bam out.bam

        samtools view -c out.bam > split_bam_umi_count.txt

        bedtools bamtobed -i out.bam -split
                | sort -k1,1 -k2,2n -k3,3n -S 5G
                > "in_bam.bed"

        bedtools map
            -a in_bam.bed
            -b exon_index
            -nonamecheck -s -f 0.95 -c 7 -o distinct -delim "|"
        | bedtools map
            -a - -b gene_index
            -nonamecheck -s -f 0.95 -c 4 -o distinct -delim "|"
        | sort -k4,4 -k2,2n -k3,3n -S 5G
        | datamash
            -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8
        | assign-reads-to-genes.py gene_index
        | awk $3 == "exonic" || $3 == "intronic" {{
                split($1, arr, "|")
                printf "%s_%s_%s	%s	%s\n", arr[3], arr[4], arr[5], $2, $3
        }}
        | sort -k2,2 -k1,1 -S 5G > in_bam.txt

        cat logfile > merge_assignment.log
        cat split_bed > key.bed
        sort -m -k1,1 -k2,2 split_gene_assign > key_ga.txt

        datamash -g 1,2 count 2 < key_ga.txt
        | gzip > key.gz

    Process stats:
        remove_dups starting reads: 6884976
        remove_dups ending reads  : 6091023


        Read assignments:
            intronic                2876079
            exonic                  1463779

** End processes 'remove duplicates, assign_genes, merge_assignment' at: Thu Jun 18 06:59:06 PDT 2020

** Start process 'count_umis_by_sample' at: Thu Jun 18 06:59:13 PDT 2020

    Process versions:
        Python 3.6.4

    Process command:
        tabulate_per_cell_counts.py
            --gene_assignment_files Barnyard_ga.txt
            --all_counts_file Barnyard.UMIs.per.cell.barcode.txt
            --intron_counts_file Barnyard.UMIs.per.cell.barcode.intronic.txt

    Process stats:
        Total cells                            : 65009
        Total cells > 100 reads                : 6775
        Total cells > 1000 reads               : 989
        Total reads in cells with > 100 reads  : 4054994

** End process 'count_umis_by_sample' at: Thu Jun 18 06:59:23 PDT 2020

** Start process 'make_matrix' at: Thu Jun 18 06:59:27 PDT 2020

    Process command:
        make_matrix.py <(zcat Barnyard.gz)
            --gene_annotation "/net/bbi/vol1/data/genomes_stage/barnyard/barnyard_rna//latest.gene.annotations"
            --key "Barnyard"
        cat /net/bbi/vol1/data/genomes_stage/barnyard/barnyard_rna//latest.gene.annotations > "Barnyard.gene_annotations.txt"

** End process 'make_matrix' at: Thu Jun 18 07:00:05 PDT 2020

** Start process 'make_cds' at: Thu Jun 18 07:00:20 PDT 2020

    Process versions:
        R version 3.6.1 (2019-07-05) -- "Action of the Toes"
            monocle3 version [1] ‘0.2.1.2’

    Process command:
        make_cds.R
            "Barnyard.umi_counts.mtx"
            "Barnyard.cell_annotations.txt"
            "Barnyard.gene_annotations.txt"
            "/net/bbi/vol1/data/genomes_stage/barnyard/barnyard_rna//latest.genes.bed"
            "Barnyard"
            "100"

** End process 'make_cds' at: Thu Jun 18 07:00:40 PDT 2020

** Start process 'apply_garnett' at: Thu Jun 18 07:00:43 PDT 2020

No Garnett classifier provided for this sample

** End process 'apply_garnett' at: Thu Jun 18 07:00:53 PDT 2020

** Start process 'run_scrublet' at: Thu Jun 18 07:00:58 PDT 2020

    Process versions:
        Python 3.6.4
            scrublet  0.2.1

    Process command:
        run_scrublet.py --key Barnyard --mat Barnyard_for_scrub.mtx

** End process 'run_scrublet' at: Thu Jun 18 07:01:51 PDT 2020

** Start processes to generate qc metrics and dashboard at: Thu Jun 18 07:01:51 PDT 2020


** End processes generate qc metrics and dashboard at: Thu Jun 18 09:15:19 PDT 2020

***** END PIPELINE *****:

***** PIPELINE READ STATS *****:

             Process       Starting reads         Ending reads               % lost      % of total lost
========================================================================================================
            Trimming              9056751              8565855                 5.42                 5.42
           Alignment              8565855              7632398                10.90                10.31
           Filtering              7632398              6884976                 9.79                 8.25
       Deduplication              6884976              6091023                11.53                 8.77
     Gene assignment              6091023              4339858                28.75                19.34

Alignment details:
                                         Count              Percent
========================================================================================================
   Total reads processed:              8565855
   Reads uniquely mapped:              6884976                80.38
      Reads multi-mapped:               747422                 8.73
         Reads too short:               598672                 6.99
`,
"293T" :  `
BBI bbi-sci Pipeline Log

Nextflow version: 20.01.0
Pipeline version: 2.0.2
Git Repository, Version, Commit ID, Session ID: https://github.com/bbi-lab/bbi-sci.git, master, daf24eb76a7927c2a7e3acb9235add93485096f9, f09594ff-4462-420b-8770-96f919ad70f6

Command:
nextflow run bbi-sci -c experiment.config -resume

***** PARAMETERS *****:

    params.run_dir:               /net/bbi/vol1/seq/nextseq/200616_NS500773_0395_AHCGGYBGXF
    params.output_dir:            /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a
    params.sample_sheet:          /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/SampleSheet.csv
    params.p7_rows:               A B G H
    params.p5_cols:               8 9 7 6
    params.demux_out:             /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/demux_out
    params.level:                 3
    params.max_cores:             16
    params.samples:               false
    params.star_file:             /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/star_file.txt
    params.gene_file:             /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/gene_file.txt
    params.umi_cutoff:            100
    params.rt_barcode_file:       default
    params.hash_list:             false
    params.max_wells_per_sample:  20


Run started at: Thu Jun 18 05:31:12 PDT 2020

***** BEGIN PIPELINE *****:

** Start process 'check_sample_sheet' at: Thu Jun 18 05:31:12 PDT 2020

    Process versions:
        Python 3.6.4

    Process command:
        check_sample_sheet.py
            --sample_sheet /net/bbi/vol1/data/sciRNAseq/nextseq_runs/RNA3-012-a/SampleSheet.csv
            --star_file /net/gs/vol1/home/hpliner/.nextflow/assets/bbi-lab/bbi-sci/bin/star_file.txt
            --level 3 --rt_barcode_file default
            --max_wells_per_samp 20

** End process 'check_sample_sheet' at: Thu Jun 18 05:31:15 PDT 2020

** Start process 'trim_fastqs' for 293T-L002.fastq.gz at: Thu Jun 18 05:31:42 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore 293T-L002.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: 293T-L002.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA 293T-L002.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 488.58 s (36 us/read; 1.69 M reads/minute).

=== Summary ===

Total reads processed:              13,753,838
Reads with adapters:                 9,264,971 (67.4%)
Reads written (passing filters):    13,753,838 (100.0%)

Total basepairs processed: 1,375,383,800 bp
Quality-trimmed:              34,631,793 bp (2.5%)
Total written (filtered):  1,028,341,257 bp (74.8%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 9264971 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 46.9%
  G: 16.2%
  T: 36.8%
  none/other: 0.2%


RUN STATISTICS FOR INPUT FILE: 293T-L002.fastq.gz
=============================================
13753838 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	789531 (5.7%)

** End process 'trim_fastqs' at: Thu Jun 18 05:41:09 PDT 2020

** Start process 'trim_fastqs' for 293T-L004.fastq.gz at: Thu Jun 18 05:31:33 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore 293T-L004.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: 293T-L004.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA 293T-L004.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 508.99 s (35 us/read; 1.72 M reads/minute).

=== Summary ===

Total reads processed:              14,554,797
Reads with adapters:                 9,747,218 (67.0%)
Reads written (passing filters):    14,554,797 (100.0%)

Total basepairs processed: 1,455,479,700 bp
Quality-trimmed:              39,544,848 bp (2.7%)
Total written (filtered):  1,082,796,678 bp (74.4%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 9747218 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 47.4%
  G: 15.8%
  T: 36.7%
  none/other: 0.1%


RUN STATISTICS FOR INPUT FILE: 293T-L004.fastq.gz
=============================================
14554797 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	863774 (5.9%)

** End process 'trim_fastqs' at: Thu Jun 18 05:41:27 PDT 2020

** Start process 'trim_fastqs' for 293T-L003.fastq.gz at: Thu Jun 18 05:31:43 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore 293T-L003.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: 293T-L003.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA 293T-L003.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 501.02 s (35 us/read; 1.71 M reads/minute).

=== Summary ===

Total reads processed:              14,298,365
Reads with adapters:                 9,570,134 (66.9%)
Reads written (passing filters):    14,298,365 (100.0%)

Total basepairs processed: 1,429,836,500 bp
Quality-trimmed:              38,033,321 bp (2.7%)
Total written (filtered):  1,065,563,756 bp (74.5%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 9570134 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 47.5%
  G: 15.8%
  T: 36.6%
  none/other: 0.1%


RUN STATISTICS FOR INPUT FILE: 293T-L003.fastq.gz
=============================================
14298365 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	840436 (5.9%)

** End process 'trim_fastqs' at: Thu Jun 18 05:41:30 PDT 2020

** Start process 'trim_fastqs' for 293T-L001.fastq.gz at: Thu Jun 18 05:31:46 PDT 2020

    Process versions:
        Python 2.7.3
        trim_galore version 0.4.1
        cutadapt version 1.9.dev2

    Process command:
        trim_galore 293T-L001.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:

SUMMARISING RUN PARAMETERS
==========================
Input filename: 293T-L001.fastq.gz
Trimming mode: single-end
Trim Galore version: 0.4.1
Cutadapt version: 1.9.dev2
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AAAAAAAA' ()
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp
All Read 1 sequences will be trimmed by 1 bp from their 3' end to avoid poor qualities or biases
Output file will be GZIP compressed


This is cutadapt 1.9.dev2 with Python 2.7.3
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AAAAAAAA 293T-L001.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 496.36 s (35 us/read; 1.70 M reads/minute).

=== Summary ===

Total reads processed:              14,045,082
Reads with adapters:                 9,379,670 (66.8%)
Reads written (passing filters):    14,045,082 (100.0%)

Total basepairs processed: 1,404,508,200 bp
Quality-trimmed:              38,503,924 bp (2.7%)
Total written (filtered):  1,047,651,677 bp (74.6%)

=== Adapter 1 ===

Sequence: AAAAAAAA; Type: regular 3'; Length: 8; Trimmed: 9379670 times.

No. of allowed errors:
0-8 bp: 0

Bases preceding removed adapters:
  A: 0.0%
  C: 47.5%
  G: 15.8%
  T: 36.7%
  none/other: 0.1%


RUN STATISTICS FOR INPUT FILE: 293T-L001.fastq.gz
=============================================
14045082 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	820913 (5.8%)

** End process 'trim_fastqs' at: Thu Jun 18 05:41:28 PDT 2020

** Start process 'align_reads' for 293T-L002_trimmed.fq.gz at: Thu Jun 18 05:52:55 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/human/human_star
            --readFilesIn 293T-L002_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/293T-L002 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:09.07.59
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
        GTF download date:   2020.04.23:09.58.35


    Process output:
                                 Started job on |	Jun 18 05:52:55
                             Started mapping on |	Jun 18 05:56:59
                                    Finished on |	Jun 18 05:58:59
       Mapping speed, Million of reads per hour |	388.93

                          Number of input reads |	12964307
                      Average input read length |	77
                                    UNIQUE READS:
                   Uniquely mapped reads number |	10282549
                        Uniquely mapped reads % |	79.31%
                          Average mapped length |	79.42
                       Number of splices: Total |	119596
            Number of splices: Annotated (sjdb) |	75103
                       Number of splices: GT/AG |	108599
                       Number of splices: GC/AG |	10513
                       Number of splices: AT/AC |	168
               Number of splices: Non-canonical |	316
                      Mismatch rate per base, % |	1.11%
                         Deletion rate per base |	0.04%
                        Deletion average length |	1.50
                        Insertion rate per base |	0.03%
                       Insertion average length |	1.11
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	1136380
             % of reads mapped to multiple loci |	8.77%
        Number of reads mapped to too many loci |	193310
             % of reads mapped to too many loci |	1.49%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	7.75%
                     % of reads unmapped: other |	2.68%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:58:59 PDT 2020

** Start process 'align_reads' for 293T-L004_trimmed.fq.gz at: Thu Jun 18 05:52:57 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/human/human_star
            --readFilesIn 293T-L004_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/293T-L004 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:09.07.59
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
        GTF download date:   2020.04.23:09.58.35


    Process output:
                                 Started job on |	Jun 18 05:52:58
                             Started mapping on |	Jun 18 05:56:59
                                    Finished on |	Jun 18 05:58:59
       Mapping speed, Million of reads per hour |	410.73

                          Number of input reads |	13691023
                      Average input read length |	77
                                    UNIQUE READS:
                   Uniquely mapped reads number |	10872659
                        Uniquely mapped reads % |	79.41%
                          Average mapped length |	79.29
                       Number of splices: Total |	125648
            Number of splices: Annotated (sjdb) |	79600
                       Number of splices: GT/AG |	113945
                       Number of splices: GC/AG |	11203
                       Number of splices: AT/AC |	175
               Number of splices: Non-canonical |	325
                      Mismatch rate per base, % |	1.00%
                         Deletion rate per base |	0.04%
                        Deletion average length |	1.50
                        Insertion rate per base |	0.03%
                       Insertion average length |	1.11
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	1199142
             % of reads mapped to multiple loci |	8.76%
        Number of reads mapped to too many loci |	206876
             % of reads mapped to too many loci |	1.51%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	7.57%
                     % of reads unmapped: other |	2.75%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 05:58:59 PDT 2020

** Start process 'align_reads' for 293T-L003_trimmed.fq.gz at: Thu Jun 18 05:58:46 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/human/human_star
            --readFilesIn 293T-L003_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/293T-L003 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:09.07.59
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
        GTF download date:   2020.04.23:09.58.35


    Process output:
                                 Started job on |	Jun 18 05:58:46
                             Started mapping on |	Jun 18 06:04:41
                                    Finished on |	Jun 18 06:06:04
       Mapping speed, Million of reads per hour |	583.72

                          Number of input reads |	13457929
                      Average input read length |	77
                                    UNIQUE READS:
                   Uniquely mapped reads number |	10714297
                        Uniquely mapped reads % |	79.61%
                          Average mapped length |	79.39
                       Number of splices: Total |	125205
            Number of splices: Annotated (sjdb) |	79540
                       Number of splices: GT/AG |	113440
                       Number of splices: GC/AG |	11234
                       Number of splices: AT/AC |	198
               Number of splices: Non-canonical |	333
                      Mismatch rate per base, % |	0.98%
                         Deletion rate per base |	0.04%
                        Deletion average length |	1.50
                        Insertion rate per base |	0.03%
                       Insertion average length |	1.11
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	1170628
             % of reads mapped to multiple loci |	8.70%
        Number of reads mapped to too many loci |	204241
             % of reads mapped to too many loci |	1.52%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	7.43%
                     % of reads unmapped: other |	2.74%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 06:06:04 PDT 2020

** Start process 'align_reads' for 293T-L001_trimmed.fq.gz at: Thu Jun 18 05:57:27 PDT 2020

    Process versions:
        STAR_2.5.2b

    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/human/human_star
            --readFilesIn 293T-L001_trimmed.fq.gz --readFilesCommand zcat
            --outFileNamePrefix ./align_out/293T-L001 --outSAMtype BAM Unsorted
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
        Genome fasta file URL: ftp.ensembl.org:/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
        FASTA download date:   2020.04.23:09.07.59
        Non REF sequences removed.

        Genome GTF file URL: ftp.ensembl.org:/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
        GTF download date:   2020.04.23:09.58.35


    Process output:
                                 Started job on |	Jun 18 05:57:27
                             Started mapping on |	Jun 18 06:09:44
                                    Finished on |	Jun 18 06:11:14
       Mapping speed, Million of reads per hour |	528.97

                          Number of input reads |	13224169
                      Average input read length |	77
                                    UNIQUE READS:
                   Uniquely mapped reads number |	10531807
                        Uniquely mapped reads % |	79.64%
                          Average mapped length |	79.43
                       Number of splices: Total |	122060
            Number of splices: Annotated (sjdb) |	77560
                       Number of splices: GT/AG |	110608
                       Number of splices: GC/AG |	10935
                       Number of splices: AT/AC |	183
               Number of splices: Non-canonical |	334
                      Mismatch rate per base, % |	0.98%
                         Deletion rate per base |	0.04%
                        Deletion average length |	1.50
                        Insertion rate per base |	0.03%
                       Insertion average length |	1.11
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	1150653
             % of reads mapped to multiple loci |	8.70%
        Number of reads mapped to too many loci |	199930
             % of reads mapped to too many loci |	1.51%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	7.41%
                     % of reads unmapped: other |	2.74%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%

** End process 'align_reads' at: Thu Jun 18 06:11:15 PDT 2020

** Start process 'sort_and_filter' for 293T-L002Aligned.out.bam at: Thu Jun 18 05:59:04 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 '293T-L002Aligned.out.bam'
            | samtools sort -@ 10 - > '293T-L002.bam'

    Process stats:
        sort_and_filter starting reads: 11418929
        sort_and_filter ending reads  : 10282549

** End process 'sort_and_filter' at: Thu Jun 18 06:03:35 PDT 2020

** Start process 'sort_and_filter' for 293T-L004Aligned.out.bam at: Thu Jun 18 05:59:05 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 '293T-L004Aligned.out.bam'
            | samtools sort -@ 10 - > '293T-L004.bam'

    Process stats:
        sort_and_filter starting reads: 12071801
        sort_and_filter ending reads  : 10872659

** End process 'sort_and_filter' at: Thu Jun 18 06:03:47 PDT 2020

** Start process 'sort_and_filter' for 293T-L003Aligned.out.bam at: Thu Jun 18 06:06:08 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 '293T-L003Aligned.out.bam'
            | samtools sort -@ 10 - > '293T-L003.bam'

    Process stats:
        sort_and_filter starting reads: 11884925
        sort_and_filter ending reads  : 10714297

** End process 'sort_and_filter' at: Thu Jun 18 06:10:29 PDT 2020

** Start process 'sort_and_filter' for 293T-L001Aligned.out.bam at: Thu Jun 18 06:11:18 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools view -bh -q 30 -F 4 '293T-L001Aligned.out.bam'
            | samtools sort -@ 10 - > '293T-L001.bam'

    Process stats:
        sort_and_filter starting reads: 11682460
        sort_and_filter ending reads  : 10531807

** End process 'sort_and_filter' at: Thu Jun 18 06:15:40 PDT 2020

** Start process 'merge_bams' at: Thu Jun 18 06:21:05 PDT 2020

    Process versions:
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.

    Process command:
        samtools merge 293T.bam 293T-L002.bam 293T-L004.bam 293T-L003.bam 293T-L001.bam

** End process 'merge_bams' at: Thu Jun 18 06:23:16 PDT 2020

** Start processes 'remove duplicates, assign_genes, merge_assignment' at: Thu Jun 18 06:23:27 PDT 2020

    Process versions:
        bedtools v2.26.0
        samtools 1.4 Using htslib 1.4 Copyright (C) 2017 Genome Research Ltd.
        bamtools 2.2.3
        Python 3.6.4

    Process command:
        mkdir split_bams
        bamtools split -in 293T.bam -reference -stub split_bams/split

        rmdup.py --bam in_bam --output_bam out.bam

        samtools view -c out.bam > split_bam_umi_count.txt

        bedtools bamtobed -i out.bam -split
                | sort -k1,1 -k2,2n -k3,3n -S 5G
                > "in_bam.bed"

        bedtools map
            -a in_bam.bed
            -b exon_index
            -nonamecheck -s -f 0.95 -c 7 -o distinct -delim "|"
        | bedtools map
            -a - -b gene_index
            -nonamecheck -s -f 0.95 -c 4 -o distinct -delim "|"
        | sort -k4,4 -k2,2n -k3,3n -S 5G
        | datamash
            -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8
        | assign-reads-to-genes.py gene_index
        | awk $3 == "exonic" || $3 == "intronic" {{
                split($1, arr, "|")
                printf "%s_%s_%s	%s	%s\n", arr[3], arr[4], arr[5], $2, $3
        }}
        | sort -k2,2 -k1,1 -S 5G > in_bam.txt

        cat logfile > merge_assignment.log
        cat split_bed > key.bed
        sort -m -k1,1 -k2,2 split_gene_assign > key_ga.txt

        datamash -g 1,2 count 2 < key_ga.txt
        | gzip > key.gz

    Process stats:
        remove_dups starting reads: 42401312
        remove_dups ending reads  : 37424453


        Read assignments:
            exonic                  6389408
            intronic               19495160

** End processes 'remove duplicates, assign_genes, merge_assignment' at: Thu Jun 18 07:00:25 PDT 2020

** Start process 'count_umis_by_sample' at: Thu Jun 18 07:00:32 PDT 2020

    Process versions:
        Python 3.6.4

    Process command:
        tabulate_per_cell_counts.py
            --gene_assignment_files 293T_ga.txt
            --all_counts_file 293T.UMIs.per.cell.barcode.txt
            --intron_counts_file 293T.UMIs.per.cell.barcode.intronic.txt

    Process stats:
        Total cells                            : 311214
        Total cells > 100 reads                : 48385
        Total cells > 1000 reads               : 5057
        Total reads in cells with > 100 reads  : 24181555

** End process 'count_umis_by_sample' at: Thu Jun 18 07:01:28 PDT 2020

** Start process 'make_matrix' at: Thu Jun 18 07:01:33 PDT 2020

    Process command:
        make_matrix.py <(zcat 293T.gz)
            --gene_annotation "/net/bbi/vol1/data/genomes_stage/human/human_rna//latest.gene.annotations"
            --key "293T"
        cat /net/bbi/vol1/data/genomes_stage/human/human_rna//latest.gene.annotations > "293T.gene_annotations.txt"

** End process 'make_matrix' at: Thu Jun 18 07:04:19 PDT 2020

** Start process 'make_cds' at: Thu Jun 18 07:04:24 PDT 2020

    Process versions:
        R version 3.6.1 (2019-07-05) -- "Action of the Toes"
            monocle3 version [1] ‘0.2.1.2’

    Process command:
        make_cds.R
            "293T.umi_counts.mtx"
            "293T.cell_annotations.txt"
            "293T.gene_annotations.txt"
            "/net/bbi/vol1/data/genomes_stage/human/human_rna//latest.genes.bed"
            "293T"
            "100"

** End process 'make_cds' at: Thu Jun 18 07:05:17 PDT 2020

** Start process 'apply_garnett' at: Thu Jun 18 07:05:24 PDT 2020

No Garnett classifier provided for this sample

** End process 'apply_garnett' at: Thu Jun 18 07:05:36 PDT 2020

** Start process 'run_scrublet' at: Thu Jun 18 07:05:42 PDT 2020

    Process versions:
        Python 3.6.4
            scrublet  0.2.1

    Process command:
        run_scrublet.py --key 293T --mat 293T_for_scrub.mtx

** End process 'run_scrublet' at: Thu Jun 18 08:00:34 PDT 2020

** Start processes to generate qc metrics and dashboard at: Thu Jun 18 08:00:34 PDT 2020


** End processes generate qc metrics and dashboard at: Thu Jun 18 09:15:19 PDT 2020

***** END PIPELINE *****:

***** PIPELINE READ STATS *****:

             Process       Starting reads         Ending reads               % lost      % of total lost
========================================================================================================
            Trimming             56652082             53337428                 5.85                 5.85
           Alignment             53337428             47058115                11.77                11.08
           Filtering             47058115             42401312                 9.90                 8.22
       Deduplication             42401312             37424453                11.74                 8.78
     Gene assignment             37424453             25884568                30.84                20.37

Alignment details:
                                         Count              Percent
========================================================================================================
   Total reads processed:             53337428
   Reads uniquely mapped:             42401312                79.50
      Reads multi-mapped:              4656803                 8.73
         Reads too short:              4020979                 7.54
`
}
