version 1.0

task germline {
    input {
        #fq2bam
        File inputFASTQ_1
        File inputFASTQ_2
        File inputRefTarball
        File? inputRecal

        String? readGroup_sampleName = "SAMPLE"
        String? readGroup_libraryName = "LIB1"
        String? readGroup_ID = "RG1"
        String? readGroup_platformName = "ILLUMINA"
        String? readGroup_PU = "unit1"

        File? inputKnownSitesVCF
        File? inputKnownSitesTBI

        String tmpDir = "tmp_fq2bam"

        # haplotypecaller
        File? intervalFile
        Boolean gvcfMode = false
        Boolean useBestPractices = false
        String haplotypecallerPassthroughOptions = ""
        String annotationArgs = ""
        
        # runtimes
        File? pbLicenseBin
        String pbPATH = "pbrun"
        String pbDocker = "bio2s-images-cn-beijing.cr.volces.com/nvidia/parabricks:4.0.0-1"
        Int nGPU = 2
        String gpuModel = "Tesla-T4"
        # String gpuDriverVersion = "460.73.01"
        Int nThreads = 24
        Int gbRAM = 64
        Int diskGB = 0
        # Int runtimeMinutes = 600
        # String hpcQueue = "gpu"
        Int maxPreemptAttempts = 1
    }

    # fq2bam
    Int auto_diskGB = if diskGB == 0 then ceil(5.0 * size(inputFASTQ_1, "GB")) + ceil(5.0 * size(inputFASTQ_2, "GB")) + ceil(3.0 * size(inputRefTarball, "GB")) + ceil(size(inputKnownSitesVCF, "GB")) + 150 else diskGB
    String best_practice_args = if useBestPractices then "--bwa-options \" -Y -K 100000000 \" " else ""
    String rgID = if readGroup_sampleName == "SAMPLE" then readGroup_ID else readGroup_sampleName + "-" + readGroup_ID
    String ref = basename(inputRefTarball, ".tar")
    String outbase = basename(basename(basename(basename(inputFASTQ_1, ".gz"), ".fastq"), ".fq"), "_1")

    # haplotypecaller
    String outVCF = outbase + ".haplotypecaller" + (if gvcfMode then '.g' else '') + ".vcf"
    String quantization_band_stub = if useBestPractices then " -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 " else ""
    String quantization_qual_stub = if useBestPractices then " --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30" else ""
    String annotation_stub_base = if useBestPractices then "-G StandardAnnotation -G StandardHCAnnotation" else annotationArgs
    String annotation_stub = if useBestPractices && gvcfMode then annotation_stub_base + " -G AS_StandardAnnotation " else annotation_stub_base

    String localTarball = basename(inputRefTarball)

    command {
        set -e
        set -x
        set -o pipefail
        mkdir -p ~{tmpDir} && \
        time tar xf ~{inputRefTarball} && \
        time ~{pbPATH} fq2bam \
        --tmp-dir ~{tmpDir} \
        --in-fq ~{inputFASTQ_1} ~{inputFASTQ_2} \
        "@RG\tID:~{rgID}\tLB:~{readGroup_libraryName}\tPL:~{readGroup_platformName}\tSM:~{readGroup_sampleName}\tPU:~{readGroup_PU}" \
        ~{best_practice_args} \
        --ref ~{ref} \
        ~{"--knownSites " + inputKnownSitesVCF + " --out-recal-file " + outbase + ".pb.BQSR-REPORT.txt"} \
        --out-bam ~{outbase}.pb.bam
 
        time ~{pbPATH} haplotypecaller \
        --in-bam ~{outbase}.pb.bam \
        --ref ~{ref} \
        --out-variants ~{outVCF} \
        ~{ if defined(inputKnownSitesVCF) then "--in-recal-file " + outbase + ".pb.BQSR-REPORT.txt" else "" } \
        ~{if gvcfMode then "--gvcf " else ""} \
        #~{"--haplotypecaller-options " + haplotypecallerPassthroughOptions } \
        ~{annotation_stub} \
        ~{quantization_band_stub} \
        ~{quantization_qual_stub}
 
    }

    output {
        File outputBAM = "~{outbase}.pb.bam"
        File outputBAI = "~{outbase}.pb.bam.bai"
        File? outputBQSR = "~{outbase}.pb.BQSR-REPORT.txt"
        File outVCF = "~{outVCF}"
    }

    runtime {
        docker : "~{pbDocker}"
        disks : "local-disk ~{auto_diskGB} SSD"
        cpu : nThreads
        memory : "~{gbRAM} GB"
        hpcMemory : gbRAM
        # hpcQueue : "~{hpcQueue}"
        # hpcRuntimeMinutes : runtimeMinutes
        gpuType : "~{gpuModel}"
        gpuCount : nGPU
        # nvidiaDriverVersion : "~{gpuDriverVersion}"
        # zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}

workflow ClaraParabricks {
    call germline

    output {
        File outputBAM = germline.outputBAM
        File outputBAI = germline.outputBAI
        File? outputBQSR = germline.outputBQSR
        File outVCF = germline.outVCF
    }
}
