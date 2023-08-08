task bwa_align {
    input {
        File fastq1
        File fastq2

        Array[File] genome_indexes

        String sample_id

        Int cpu = 64
        Int memory = 128
        Int disk_size_gb = ceil(6 * ceil(size(fastq1, 'GB')) /10) * 10 + 20
        String docker_img
    }
    String reads_group = "'@RG\\tID:~{sample_id}\\tSM:~{sample_id}\\tPL:ILLUMINA'"
    String sorted_bam = sample_id + '.sorted.bam'
    String sorted_bam_index = sorted_bam + '.bai'

    command <<<
       cpu_cores=$(nproc)

       bwa mem \
            -Y -K 100000000 \
            -t $cpu_cores \
            -R ~{reads_group} \
            ~{genome_indexes[0]} \
            ~{fastq1} ~{fastq2} \
       | samtools sort -@ $cpu_cores -o ~{sorted_bam}

       samtools index ~{sorted_bam}

    >>>

    runtime {
        cpu: "~{cpu}"
        memory: '~{memory} GB'
        disk: "~{disk_size_gb} GB"
        docker: docker_img
    }

    output {
        Pair[File,File] sorted_bam_output = (sorted_bam, sorted_bam_index)
    }
}

workflow wgs {

    input {
        String sample_id
    
        File fastq1
        File fastq2

        Boolean if_down_fastq
        Int? target_depth
        Float? origin_depth

        Boolean phred64 
        Boolean fix_mgi_id
        String? adapter_sequence
        String? adapter_sequence_r2
        Int? reads_to_process 

        Array[File] genome_indexes

        String docker_img_germline_tools = 'cr-cn-beijing.ivolces.com/beeding/germline_tools:v3'
    }
    
    call tasks.fastp_clean {
        input:
            sample_id = sample_id,
            fastq1 = run_fastq1,
            fastq2 = run_fastq2,
            phred64 = phred64,
            fix_mgi_id = fix_mgi_id,
            adapter_sequence = adapter_sequence,
            adapter_sequence_r2 = adapter_sequence_r2,
            reads_to_process = reads_to_process,
            docker_img = docker_img_germline_tools
    }