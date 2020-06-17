version 1.0

workflow strelka {
  input {
    File tumorBam
    File tumorBai
    File normalBam
    File normalBai
    String outputFileNamePrefix
    Int isSkipDepthFilters
    Int maxInputDepth
    Boolean doBamSort
  }

  parameter_meta {
    tumorBam: "Tumor bam file"
    normalBam: "Matched normal bam file"
    isSkipDepthFilters: "Used to generated config.ini. Value of 1 to skip depth filtration for whole exome or other targeted sequencing data"
    maxInputDepth: "Used to generated config.ini. Strelka will not accept input reads above this depth."
    doBamSort: "Flag to indicate whether bam should be sorted. Default for Strelka: true."
  }

  call generateConfig {
    input:
      isSkipDepthFilters = isSkipDepthFilters,
      maxInputDepth = maxInputDepth
  }

  if (doBamSort) {
    call sortBams {
      input:
        tumorBam = tumorBam,
        normalBam = normalBam
    }
  }

  call runStrelka {
    input:
      tumorBam = select_first([sortBams.tumorBamSorted, tumorBam]),
      tumorBai = select_first([sortBams.tumorBaiReordered, tumorBai]),
      normalBam = select_first([sortBams.normalBamSorted, normalBam]),
      normalBai = select_first([sortBams.normalBaiReordered, normalBai]),
      outputFileNamePrefix = outputFileNamePrefix,
      configIni = generateConfig.configIni
  }

  String normal = basename(normalBam, '.bam')
  String tumor = basename(tumorBam, '.bam')
  String header = "NORMAL:" + normal + " TUMOR:" + tumor

  call updateVcfHeader {
    input:
      snv = runStrelka.snv,
      indel = runStrelka.indel,
      header = header
  }

  output {
    File snv = updateVcfHeader.updatedSnv
    File snvIndex = updateVcfHeader.snvIndex
    File indel = updateVcfHeader.updatedIndel
    File indelIndex = updateVcfHeader.indelIndex
  }

  meta {
    author: "Angie Mosquera"
    email: "angie.mosquera@oicr.on.ca"
    description: "Workflow for Strelka somatic small-variant caller."
    dependencies: [
      {
        name: "samtools/1.9",
        url: "http://www.htslib.org/"
      },
      {
        name: "picard/2.21.2",
        url: "https://github.com/broadinstitute/picard/releases"
      },
      {
        name: "strelka/1.0.15",
        url: "https://sites.google.com/site/strelkasomaticvariantcaller/home/download"
      },
      {
        name: "tabix/1.9",
        url: "http://www.htslib.org/"
      }
    ]
  }
}

task generateConfig {
  input {
    String modules = ""
    Int isSkipDepthFilters
    Int maxInputDepth
    Float depthFilterMultiple = 3.0
    Float snvMaxFilteredBasecallFrac = 0.4
    Float snvMaxSpanningDeletionFrac = 0.75
    Int indelMaxRefRepeat = 8
    Float indelMaxWindowFilteredBasecallFrac = 0.3
    Int indelMaxIntHpolLength = 14
    Float ssnvPrior = 0.000001
    Float sindelPrior = 0.000001
    Float ssnvNoise = 0.0000005
    Float sindelNoise = 0.0000001
    Float ssnvNoiseStrandBiasFrac = 0.5
    Int minTier1Mapq = 20
    Int minTier2Mapq = 5
    Int ssnvQuality_LowerBound = 15
    Int sindelQuality_LowerBound = 30
    Int isWriteRealignedBam = 0
    Int binSize = 25000000
    String? extraStrelkaArguments
    Int memory = 4
    Int timeout = 1
  }

  parameter_meta {
    modules: "Names and versions of modules to load. None in this case."
    isSkipDepthFilters: "Flag to skip depth filtration. 0 for WG, 1 for EX or TS."
    maxInputDepth: "Strelka will not accept input reads above this depth."
    depthFilterMultiple: "If depth filter is not skipped, all variants which occur at a depth greather than depthFilterMultiple*chromosome mean depth will be filtered out."
    snvMaxFilteredBasecallFrac: "Somatic SNV calls are filtered at sites where greater than this fraction of basecalls have been removed by the mismatch density filter in either sample."
    snvMaxSpanningDeletionFrac: "Somatic SNV calls are filtered at sites where greater than this fraction of overlapping reads contain deletions which span the SNV call site."
    indelMaxRefRepeat: "Somatic indel calls are filtered if they represent an expansion or contraction of a repeated pattern with a repeat count greater than indelMaxRefRepeat in the reference"
    indelMaxWindowFilteredBasecallFrac: "Somatic indel calls are filtered if greater than this fraction of basecalls in a window extending 50 bases to each side of an indel's call position have been removed by the mismatch density filter."
    indelMaxIntHpolLength: "Somatic indels are filtered if they overlap ’interrupted homopolymers’ greater than this length. The term 'interrupted homopolymer' is used to indicate the longest homopolymer which can be found intersecting or adjacent to the called indel when a single non-homopolymer base is allowed."
    ssnvPrior: "Prior probability of a somatic snv."
    sindelPrior: "Prior probability of a somatic indel."
    ssnvNoise: "Probability of an snv noise allele."
    sindelNoise: "Probability of an indel noise allele."
    ssnvNoiseStrandBiasFrac: "Fraction of snv noise attributed to strand-bias."
    minTier1Mapq: "Minimum MAPQ score for PE reads at tier1."
    minTier2Mapq: "Minimum MAPQ score for PE and SE reads at tier2."
    ssnvQuality_LowerBound: "Somatic quality score (QSS_NT, NT=ref) below which somatic SNVs are marked as filtered."
    sindelQuality_LowerBound: "Somatic quality score (QSI_NT, NT=ref) below which somatic indels are marked as filtered"
    isWriteRealignedBam: "Optionally write out read alignments which were altered during the realignment step"
    binSize: "Jobs are parallelized over segments of the reference genome no larger than this size"
    extraStrelkaArguments: "Extra arguments for Strelka."
    memory: "Memory allocated for this job."
    timeout: "Hours before task timeout."
  }

  command <<<
  set -euo pipefail

  touch config.ini

  configText="; User configuration options for Strelka somatic small-variant caller workflow:

  [user]

  isSkipDepthFilters = ~{isSkipDepthFilters}

  maxInputDepth = ~{maxInputDepth}

  depthFilterMultiple = ~{depthFilterMultiple}

  snvMaxFilteredBasecallFrac = ~{snvMaxFilteredBasecallFrac}

  snvMaxSpanningDeletionFrac = ~{snvMaxSpanningDeletionFrac}

  indelMaxRefRepeat = ~{indelMaxRefRepeat}

  indelMaxWindowFilteredBasecallFrac = ~{indelMaxWindowFilteredBasecallFrac}

  indelMaxIntHpolLength = ~{indelMaxIntHpolLength}

  ssnvPrior = ~{ssnvPrior}

  sindelPrior = ~{sindelPrior}

  ssnvNoise = ~{ssnvNoise}

  sindelNoise = ~{sindelNoise}

  ssnvNoiseStrandBiasFrac = ~{ssnvNoiseStrandBiasFrac}

  minTier1Mapq = ~{minTier1Mapq}

  minTier2Mapq = ~{minTier2Mapq}

  ssnvQuality_LowerBound = ~{ssnvQuality_LowerBound}

  sindelQuality_LowerBound = ~{sindelQuality_LowerBound}

  isWriteRealignedBam = ~{isWriteRealignedBam}

  binSize = ~{binSize}"

  if [ -z "~{extraStrelkaArguments}" ]; then
    configText+="; Additional arguments passed to strelka.

     extraStrelkaArguments = ~{extraStrelkaArguments}"
  fi

  echo "$configText" >> config.ini
  >>>

  runtime {
    modules: "~{modules}"
    memory: "~{memory} GB"
    timeout: "~{timeout}"
  }

  output {
    File configIni = "config.ini"
  }

  meta {
    output_meta: {
      configIni: "config.ini file to be used in Strelka call."
    }
  }
}

task sortBams {
  input {
    String modules = "samtools/1.9 hg19/p13 picard/2.21.2"
    File tumorBam
    File normalBam
    Int picardMaxMemMb = 10000
    String refDict = "$HG19_ROOT/hg19_random.dict"
    String picard = "$PICARD_ROOT/picard.jar"
    Int memory = 16
    Int timeout = 4
    Int threads = 8
  }

  parameter_meta {
    modules: "Names and versions of modules to load."
    tumorBam: "Tumor bam file to be sorted."
    normalBam: "Normal bam file to be sorted."
    picardMaxMemMb: "Max amount of memory to be used by Picard ReorderSam"
    refDict: "Reference sequence dictionary file."
    memory: "Memory allocated for this job."
    timeout: "Hours before task timeout."
  }

  String tumorBamBasename = basename(tumorBam, '.bam')
  String normalBamBasename = basename(normalBam, '.bam')

  command <<<
    set -euo pipefail

    samtools sort ~{tumorBam} -o tumor.sort.bam
    samtools sort ~{normalBam} -o normal.sort.bam

    java -Xmx~{picardMaxMemMb}M \
    -jar ~{picard} ReorderSam \
    INPUT=tumor.sort.bam \
    OUTPUT="~{tumorBamBasename}.sort.reordered.bam" \
    SEQUENCE_DICTIONARY=~{refDict}

    java -Xmx~{picardMaxMemMb}M \
    -jar ~{picard} ReorderSam \
    INPUT=normal.sort.bam \
    OUTPUT="~{normalBamBasename}.sort.reordered.bam" \
    SEQUENCE_DICTIONARY=~{refDict}

    samtools index "~{tumorBamBasename}.sort.reordered.bam" "~{tumorBamBasename}.sort.reordered.bam.bai"
    samtools index "~{normalBamBasename}.sort.reordered.bam" "~{normalBamBasename}.sort.reordered.bam.bai"
  >>>

  runtime {
    modules: "~{modules}"
    memory: "~{memory} GB"
    timeout: "~{timeout}"
    cpu: "~{threads}"
  }

  output {
    File tumorBamSorted = "~{tumorBamBasename}.sort.reordered.bam"
    File tumorBaiReordered = "~{tumorBamBasename}.sort.reordered.bam.bai"
    File normalBamSorted = "~{normalBamBasename}.sort.reordered.bam"
    File normalBaiReordered = "~{normalBamBasename}.sort.reordered.bam.bai"
  }

  meta {
    output_meta: {
      tumorBamSorted: "Sorted tumor bam file.",
      normalBamSorted: "Sorted normal bam file."
    }
  }
}

task runStrelka {
  input {
    String modules = "strelka/1.0.15 hg19/p13"
    File tumorBam
    File tumorBai
    File normalBam
    File normalBai
    File configIni
    String outputFileNamePrefix
    String strelkaTag = "strelka"
    String refFasta = "$HG19_ROOT/hg19_random.fa"
    Int memory = 32
    Int threads = 16
    Int timeout = 72
  }

  parameter_meta {
    modules: "Names and versions of modules to load."
    tumorBam: "Sorted or unsorted tumor bam file."
    normalBam: "Sorted or unsorted normal bam file."
    configIni: "Config file with all Streka configuration parameters."
    refFasta: "Reference FASTA file."
    memory: "Memory allocated for this job."
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout."
  }

  String outputSnv = outputFileNamePrefix + ".final." + strelkaTag + ".snvs.vcf"
  String outputIndel = outputFileNamePrefix + ".final." + strelkaTag + ".indels.vcf"

  command <<<
    set -euo pipefail

    dataDir=$(pwd)

    ${STRELKA_ROOT}/bin/configureStrelkaWorkflow.pl \
    --normal=~{normalBam} \
    --tumor=~{tumorBam} \
    --ref=~{refFasta} \
    --config=~{configIni} \
    --output-dir="strelkaAnalysis"

    make -C ./strelkaAnalysis

    mv $dataDir/strelkaAnalysis/results/all.somatic.snvs.vcf $dataDir/~{outputSnv}
    mv $dataDir/strelkaAnalysis/results/all.somatic.indels.vcf $dataDir/~{outputIndel}
  >>>

  runtime {
    cpu: "~{threads}"
    modules: "~{modules}"
    memory: "~{memory} GB"
    timeout: "~{timeout}"
  }

  output {
    File snv = "~{outputSnv}"
    File indel = "~{outputIndel}"
  }

  meta {
    output_meta: {
      snv: "All somatic SNV predictions.",
      indel: "All somatic indel predictions."
    }
  }
}

task updateVcfHeader {
  input {
    String modules = "update-vcf-header-deps/0.0.1 tabix/1.9"
    File snv
    File indel
    String header
    String caller = "strelka"
    String strelkaVersion = "1.0.15"
    String reference = "hg19"
    Int memory = 8
    Int timeout = 12
    Int threads = 4
  }

  String snvFileName = basename(snv)
  String indelFileName = basename(indel)

  parameter_meta {
    modules: "Names and versions of modules to load."
    snv: ""
    indel: ""
    memory: "Memory allocated for this job."
    timeout: "Hours before task timeout."
  }

  command <<<
    set -euo pipefail

    python3 <<CODE
    import vcf
    import re

    format_lines = [
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype, contructed from SGT INFO via external modification\">\n",
      "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order "
      "listed\">\n"]

    snv_reader = vcf.Reader(filename="~{snv}", compressed=False)
    indel_reader = vcf.Reader(filename="~{indel}", compressed=False)


    def generate_AD(samples, ref, alt):
      ref_ref_count = ','.join([str(samples[0][x + 'U'][0]) for x in ref])
      alt_ref_count = ','.join([str(samples[1][x + 'U'][0]) for x in ref])

      if alt[0] is None:
        ref_alt_count = "0"
        alt_alt_count = "0"
      else:
        ref_alt_count = ','.join([str(samples[0][str(x) + 'U'][0]) for x in alt])
        alt_alt_count = ','.join([str(samples[1][str(x) + 'U'][0]) for x in alt])

      ref_count = [ref_ref_count, ref_alt_count]
      alt_count = [alt_ref_count, alt_alt_count]

      return [ref_count, alt_count] if samples[0].sample == 'NORMAL' else [alt_count, ref_count]


    def modify_header_and_records(reader, is_snv):
      header_lines = ["##fileformat=" + reader.metadata['fileformat'] + "\n", "##source=~{caller}\n",
                      "##source_version=~{strelkaVersion}\n", "##reference=" + reader.metadata['reference'] + "\n"]

      for c in reader.contigs:
        if re.search('_', reader.contigs[c].id):
          continue
        header_lines.append("##contig=<ID=" + reader.contigs[c].id + ",length=" + str(reader.contigs[c].length) +
                            ",assembly=~{reference}>\n")

      for i in reader.infos:
        num = reader.infos[i].num
        if reader.infos[i].num is None:
          num = '.'
        header_lines.append("##INFO=<ID=" + reader.infos[i].id +
                            ",Number=" + str(num) + ",Type=" + reader.infos[i].type +
                            ",Description=\"" + reader.infos[i].desc + "\">\n")
      add_gt = True
      add_ad = True
      for f in reader.formats:
        num = reader.formats[f].num if reader.formats[f].num is not None else '.'
        if reader.formats[f].id == 'GT':
          add_gt = False
        if is_snv and reader.formats[f].id == 'AD':
          add_ad = False
        header_lines.append("##FORMAT=<ID=" + reader.formats[f].id +
                            ",Number=" + str(num) + ",Type=" + reader.formats[f].type +
                            ",Description=\"" + reader.formats[f].desc + "\">\n")

      if add_gt:
        header_lines.append(format_lines[0])
      if is_snv and add_ad:
        header_lines.append(format_lines[1])

      for f in reader.filters:
        header_lines.append("##FILTER=<ID=" + reader.filters[f].id +
                            ",Description=\"" + reader.filters[f].desc + "\">\n")

      if 'inputs' not in reader.metadata.keys():
        header_lines.append("##inputs=~{header}\n")

      header_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR']
      header_lines.append("#" + "\t".join(header_fields) + "\n")

      record_lines = []
      for record in reader:
        if re.search('_', record.CHROM):
          continue

        record_data = [record.CHROM, str(record.POS), '.' if record.ID is None else record.ID, record.REF,
                       ",".join(map(str, '.' if record.ALT[0] is None else record.ALT)),
                       '.' if record.QUAL is None else record.QUAL, "PASS" if len(record.FILTER) == 0 else ";".join(record.FILTER)]
        info_data = []
        for field in record.INFO:
            if isinstance(record.INFO[field], list):
                info_data.append("=".join([field, ",".join(map(str, record.INFO[field]))]))
            else:
                info_data.append("=".join([field, "." if record.INFO[field] is None else str(record.INFO[field])]))
        info_string = ";".join(info_data) if len(info_data) > 0 else "."
        record_data.append(info_string)

        format_data = []
        if add_gt:
          add_format = 'GT:'
          if add_ad and is_snv:
            add_format += 'AD'
          format_data.append(add_format)

        for field in record.FORMAT.split(":"):
          format_data.append(field)
        record_data.append(":".join(format_data))

        if add_gt:
          gt = ['0', '0/1'] if record.samples[0].sample == 'NORMAL' else ['0/1', '0']
          if add_ad and is_snv:
            ad = generate_AD(record.samples, [record.REF], record.ALT)

        sample_data = []
        for sample in record.samples:
          if add_gt:
            format_data = [gt.pop(0), ",".join(ad.pop(0))] if add_ad and is_snv else [gt.pop(0)]
          else:
            format_data = []
          for field in record.FORMAT.split(":"):
            if isinstance(sample[field], list):
              format_data.append(",".join(map(str, sample[field])))
            else:
              format_data.append("." if sample[field] is None else str(sample[field]))
          sample_data.append(":".join(format_data))

        for s in sample_data:
          record_data.append(s)
        record_lines.append("\t".join(record_data) + "\n")

      return header_lines, record_lines


    modifiedSnv = modify_header_and_records(snv_reader, True)
    modifiedIndel = modify_header_and_records(indel_reader, False)

    with open("~{snvFileName}", mode='w+') as out:
      out.writelines(modifiedSnv[0])
      out.writelines(modifiedSnv[1])
    out.close()

    with open("~{indelFileName}", mode='w+') as out:
      out.writelines(modifiedIndel[0])
      out.writelines(modifiedIndel[1])
    out.close()
    CODE

    bgzip ~{snvFileName} && tabix -p vcf ~{snvFileName}.gz
    bgzip ~{indelFileName} && tabix -p vcf ~{indelFileName}.gz
  >>>

  runtime {
    modules: "~{modules}"
    memory: "~{memory} GB"
    timeout: "~{timeout}"
    cpu: "~{threads}"
  }

  output {
    File updatedSnv = "~{snvFileName}.gz"
    File snvIndex = "~{snvFileName}.gz.tbi"
    File updatedIndel = "~{indelFileName}.gz"
    File indelIndex = "~{indelFileName}.gz.tbi"
  }

  meta {
    output_meta: {
      updatedSnv: "snv vcf with added GT and AD fields and updated header line.",
      snvIndex: "Index of bgzipped snv vcf",
      updatedIndel: "",
      indelIndex: "Index of bgzipped indel vcf"
    }
  }
}