# downsample reads
def downsample_reads(R1_file, R2_file, R1_out, R2_out, genome_size):
    
    R1_file = R1_file.pop()
    R2_file = R2_file.pop()
    R1_out = R1_out.pop()
    R2_out = R2_out.pop()

    gsize = genome_size.pop()

    print("Using Genome Size: %s to calculate coverage" % gsize)

    # Extract basic fastq reads stats with seqtk
    seqtk_check = "/nfs/esnitkin/bin_group/seqtk/seqtk fqchk -q3 %s > %s_fastqchk.txt" % (R1_file, R1_file)

    print(seqtk_check)
    
    try:
        os.system(seqtk_check)
    except sp.CalledProcessError: 
        print('Error running seqtk for extracting fastq statistics.')
        sys.exit(1)

    with open("%s_fastqchk.txt" % R1_file, 'r') as file_open: 
        for line in file_open:
            if line.startswith('min_len'):
                line_split = line.split(';')
                min_len = line_split[0].split(': ')[1]
                max_len = line_split[1].split(': ')[1]
                avg_len = line_split[2].split(': ')[1]
            if line.startswith('ALL'):
                line_split = line.split('\t')
                total_bases = int(line_split[1]) * 2
    file_open.close()

    print('Average Read Length: %s' % avg_len)

    print('Total number of bases in fastq: %s' % total_bases)

    # Calculate original depth and check if it needs to be downsampled to a default coverage.
    ori_coverage_depth = int(total_bases / gsize)
    
    print('Original Covarage Depth: %s x' % ori_coverage_depth)


    if ori_coverage_depth > 100:
        # Downsample to 100
        factor = float(100 / float(ori_coverage_depth))
        #r1_sub = "/tmp/%s" % os.path.basename(R1_file)
        r1_sub = R1_out

        # Downsample using seqtk
        try:
            print("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (R1_file, factor, r1_sub)) 
            seqtk_downsample = "/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (R1_file, factor, r1_sub) 
            os.system(seqtk_downsample)
        except sp.CalledProcessError:
            print('Error running seqtk for downsampling raw fastq reads.')
            sys.exit(1)

        if R2_file:
            r2_sub = R2_out
            try:
                print("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (R2_file, factor, r2_sub))  
                os.system("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (R2_file, factor, r2_sub))  
            except sp.CalledProcessError:
                print('Error running seqtk for downsampling raw fastq reads.')
                sys.exit(1)
        else:
            r2_sub = "None"
    else:
        r1_sub = R1_file
        r2_sub = R2_file
        os.system("cp %s %s" % (R1_file, R1_out))
        os.system("cp %s %s" % (R2_file, R2_out))

rule downsample_clean_reads:
    input:
        r1=trimmed_r1_for_downsample,
        r2=trimmed_r2_for_downsample
    output:
        outr1="results/{prefix}/downsample/{sample}/{sample}_R1_trim_paired.fastq.gz",
        outr2="results/{prefix}/downsample/{sample}/{sample}_R2_trim_paired.fastq.gz"
    params:
        gsize=config["genome_size"]
    threads: 3
    resources:
        mem_mb=1000,
        runtime=15
    benchmark: 
        "benchmarks/{prefix}/downsample_clean_reads/{sample}.benchmark.tsv"
    run:
        downsample_reads([input.r1], [input.r2], [output.outr1], [output.outr2], [params.gsize])
