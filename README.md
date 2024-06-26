# strelka

Workflow for Strelka somatic small-variant caller.

## Overview

## Dependencies

* [samtools 1.9](http://www.htslib.org/)
* [picard 2.21.2](https://github.com/broadinstitute/picard/releases)
* [strelka 1.0.15](https://sites.google.com/site/strelkasomaticvariantcaller/home/download)
* [tabix 1.9](http://www.htslib.org/)


## Usage

### Cromwell
```
java -jar cromwell.jar run strelka.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`tumorBam`|File|Tumor bam file
`tumorBai`|File|index of Tumor bam file
`normalBam`|File|Matched normal bam file
`normalBai`|File|index of Matched normal bam file
`outputFileNamePrefix`|String|File output prefix
`isSkipDepthFilters`|Int|Used to generated config.ini. Value of 1 to skip depth filtration for whole exome or other targeted sequencing data
`maxInputDepth`|Int|Used to generated config.ini. Strelka will not accept input reads above this depth.
`doBamSort`|Boolean|Flag to indicate whether bam should be sorted. Default for Strelka: true.


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`generateConfig.modules`|String|""|Names and versions of modules to load. None in this case.
`generateConfig.depthFilterMultiple`|Float|3.0|If depth filter is not skipped, all variants which occur at a depth greather than depthFilterMultiple*chromosome mean depth will be filtered out.
`generateConfig.snvMaxFilteredBasecallFrac`|Float|0.4|Somatic SNV calls are filtered at sites where greater than this fraction of basecalls have been removed by the mismatch density filter in either sample.
`generateConfig.snvMaxSpanningDeletionFrac`|Float|0.75|Somatic SNV calls are filtered at sites where greater than this fraction of overlapping reads contain deletions which span the SNV call site.
`generateConfig.indelMaxRefRepeat`|Int|8|Somatic indel calls are filtered if they represent an expansion or contraction of a repeated pattern with a repeat count greater than indelMaxRefRepeat in the reference
`generateConfig.indelMaxWindowFilteredBasecallFrac`|Float|0.3|Somatic indel calls are filtered if greater than this fraction of basecalls in a window extending 50 bases to each side of an indel's call position have been removed by the mismatch density filter.
`generateConfig.indelMaxIntHpolLength`|Int|14|Somatic indels are filtered if they overlap ’interrupted homopolymers’ greater than this length. The term 'interrupted homopolymer' is used to indicate the longest homopolymer which can be found intersecting or adjacent to the called indel when a single non-homopolymer base is allowed.
`generateConfig.ssnvPrior`|Float|1e-06|Prior probability of a somatic snv.
`generateConfig.sindelPrior`|Float|1e-06|Prior probability of a somatic indel.
`generateConfig.ssnvNoise`|Float|5e-07|Probability of an snv noise allele.
`generateConfig.sindelNoise`|Float|1e-07|Probability of an indel noise allele.
`generateConfig.ssnvNoiseStrandBiasFrac`|Float|0.5|Fraction of snv noise attributed to strand-bias.
`generateConfig.minTier1Mapq`|Int|20|Minimum MAPQ score for PE reads at tier1.
`generateConfig.minTier2Mapq`|Int|5|Minimum MAPQ score for PE and SE reads at tier2.
`generateConfig.ssnvQuality_LowerBound`|Int|15|Somatic quality score (QSS_NT, NT=ref) below which somatic SNVs are marked as filtered.
`generateConfig.sindelQuality_LowerBound`|Int|30|Somatic quality score (QSI_NT, NT=ref) below which somatic indels are marked as filtered
`generateConfig.isWriteRealignedBam`|Int|0|Optionally write out read alignments which were altered during the realignment step
`generateConfig.binSize`|Int|25000000|Jobs are parallelized over segments of the reference genome no larger than this size
`generateConfig.extraStrelkaArguments`|String?|None|Extra arguments for Strelka.
`generateConfig.memory`|Int|4|Memory allocated for this job.
`generateConfig.timeout`|Int|1|Hours before task timeout.
`sortBams.modules`|String|"samtools/1.9 hg19/p13 picard/2.21.2"|Names and versions of modules to load.
`sortBams.picardMaxMemMb`|Int|10000|Max amount of memory to be used by Picard ReorderSam
`sortBams.refDict`|String|"$HG19_ROOT/hg19_random.dict"|Reference sequence dictionary file.
`sortBams.picard`|String|"$PICARD_ROOT/picard.jar"|Path to picard jar
`sortBams.memory`|Int|16|Memory allocated for this job.
`sortBams.timeout`|Int|4|Hours before task timeout.
`sortBams.threads`|Int|8|Threads used by this task
`runStrelka.modules`|String|"strelka/1.0.15 hg19/p13"|Names and versions of modules to load.
`runStrelka.strelkaTag`|String|"strelka"|Tag used by strelka for its outputs
`runStrelka.refFasta`|String|"$HG19_ROOT/hg19_random.fa"|Reference FASTA file.
`runStrelka.memory`|Int|32|Memory allocated for this job.
`runStrelka.threads`|Int|16|Requested CPU threads
`runStrelka.timeout`|Int|72|Hours before task timeout.
`updateVcfHeader.modules`|String|"update-vcf-header-deps/0.0.1 tabix/1.9"|Names and versions of modules to load.
`updateVcfHeader.caller`|String|"strelka"|strelka
`updateVcfHeader.strelkaVersion`|String|"1.0.15"|Version of strekla running the analysis
`updateVcfHeader.reference`|String|"hg19"|Reference id, such as hg19
`updateVcfHeader.memory`|Int|8|Memory allocated for this job.
`updateVcfHeader.timeout`|Int|12|Hours before task timeout.
`updateVcfHeader.threads`|Int|4|Threads used by this task


### Outputs

Output | Type | Description | Labels
---|---|---|---
`snv`|File|SNV calls in vcf format|vidarr_label: snv
`snvIndex`|File|Index file for SNV calls in vcf format|vidarr_label: snvIndex
`indel`|File|Indel calls in vcf format|vidarr_label: indel
`indelIndex`|File|Index file for Indel calls in vcf format|vidarr_label: indelIndex


## Commands
This section lists command(s) run by strelka workflow
 
* Running strelka
 
Strelka is an older version of strelka2 (aka strelka somatic), SNV and indel caller
 
### Configure
 
```
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
   if [ ! -z "~{extraStrelkaArguments}" ]; then
     configText+="extraStrelkaArguments = ~{extraStrelkaArguments}"
   fi
   echo "$configText" >> config.ini
```
 
### Pre-process inputs
 
```
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
```
 
### Pre-flight checks and configuration, run strelka
 
```
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
```
 
### Post-process calls
 
```
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
```
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
