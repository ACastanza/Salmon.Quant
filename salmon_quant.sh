#! /bin/bash
# Using getopt
set -e

trap abort ERR PROF
abort()
{
rm -rf salmon_index

    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

# infile="/Users/acastanza/salmon_analysis/input.file.list.txt"
# index="/Users/acastanza/Downloads/gencode.vM25.annotation.salmon_index.kmer.17"

while getopts ":r:i:l:e:n:s:g:p:t:o:h:d:q:m:f:c:b:a:j:" opt; do
    case $opt in
        r)
            infile=`realpath $OPTARG`
            ;;
        i)
            index=`realpath $OPTARG`
            ;;
        l)
            lib="$OPTARG"
            ;;
        e)
            resampling="$OPTARG"
            ;;
        n)
            bootstraps="$OPTARG"
            ;;
        s)
            seqBias="$OPTARG"
            ;;
        g)
            gcBias="$OPTARG"
            ;;
        p)
            posBias="$OPTARG"
            ;;
        t)
            BOWTIE="$OPTARG"
            ;;
        o)
            recoverOrphans="$OPTARG"
            ;;
        h)
            hardFilter="$OPTARG"
            ;;
        d)
            allowDovetail="$OPTARG"
            ;;
        q)
            dumpEq="$OPTARG"
            ;;
        m)
            reduceGCMemory="$OPTARG"
            ;;
        f)
            rangeFactor="$OPTARG"
            ;;
        c)
            factorbins="$OPTARG"
            ;;
        b)
            biasSamp="$OPTARG"
            ;;
        a)
            bias="$OPTARG"
            ;;
        j)
            threads="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            abort
            ;;
    esac
done

 mkdir -p salmon_index
 tar -zxvf $index -C salmon_index

# Process Input files

 Rall=($(grep -E '.fastq.gz|.fastq.gz|.fq.gz|.fq.gz' $infile | sort))
 R1=($(grep -E '_R1.fastq.gz|_1.fastq.gz|_R1.fq.gz|_1.fq.gz' $infile | sort))
 R2=($(grep -E '_R2.fastq.gz|_2.fastq.gz|_R2.fq.gz|_2.fq.gz' $infile | sort))
 Rpaired=("${R1[@]}" "${R2[@]}")

for sample in "${Rpaired[@]}"; do
  for i in "${!Rall[@]}"; do
    if [[ ${Rall[i]} = $sample ]]; then
      unset 'Rall[i]'
    fi
  done
done
RU=()
for i in "${!Rall[@]}"; do
    RU+=( "${Rall[i]}" )
done

# Begin Sample Quantification

if [[ "${#R1[@]}" > 0 ]]; then
 if [[ "${#R1[@]}" -eq "${#R2[@]}" ]]; then
 	for (( i=0; i<"${#R1[@]}"; i++ )); do

 outdir=$(basename ${R1[$i]})
 outdir=${outdir/%_R1.fastq.gz}
 outdir=${outdir/%_1.fastq.gz}
 outdir=${outdir/%_R1.fq.gz}
 outdir=${outdir/%_1.fq.gz}

    echo -e "\n" ;
    echo -e "--Input file(s) are:\t""${R1[$i]}"",""${R2[$i]}" ;
    echo -e "--Processing in Paired-end mode" ;
    echo -e "--Output directory is:\t""${outdir}" ;
    mkdir -p "$outdir" ;


params=()
[[ $resampling == "BOOT" ]] && params+=(--numBootstraps=$bootstraps)
[[ $resampling == "GIBB" ]] && params+=(--numGibbsSamples=$bootstraps)
[[ $validateMappings == true ]] && params+=(--validateMappings)
[[ $seqBias == true ]] && params+=(--seqBias)
[[ $gcBias == true ]] && params+=(--gcBias)
[[ $posBias == true ]] && params+=(--posBias)
[[ $BOWTIE == "BT2" ]] && params+=(--mimicBT2)
[[ $BOWTIE == "Strict" ]] && params+=(--mimicStrictBT2)
[[ $recoverOrphans == true ]] && params+=(--recoverOrphans)
[[ $hardFilter == true ]] && params+=(--hardFilter)
[[ $allowDovetail == true ]] && params+=(--allowDovetail)
[[ $dumpEq == true ]] && params+=(--dumpEq)
[[ $reduceGCMemory == true ]] && params+=(--reduceGCMemory)
[[ $rangeFactor == true ]] && params+=(--rangeFactorizationBins=$factorbins)
[[ $biasSamp == true ]] && params+=(--biasSpeedSamp=$bias)

# [[ $CONDITION == true ]] && params+=(--param)

salmon quant \
      --no-version-check \
      --index="salmon_index" \
      --threads=$threads \
      --libType=$lib \
      -1 "${R1[$i]}" \
      -2 "${R2[$i]}" \
      "${params[@]}" \
      --output=$outdir ;

cp $outdir/quant.sf $outdir.quant.sf
tar -czvf $outdir.salmon_quant.tar.gz -C $outdir .
rm -rf $outdir

    echo $outdir": Done." ;
  done
 else
 echo "Paired-end files detected but count of Read 1 ("${#R1[@]}") didn't match count of Read 2 ("${#R2[@]}")"
 echo ""${#RU[@]}" possible single-end file(s) detected"
 fi
fi

if [[ "${#RU[@]}" > 0 ]]; then
 for (( i=0; i<"${#RU[@]}"; i++ )); do

 outdir=$(basename ${R1[$i]})
 outdir=${outdir/%.fastq.gz}
 outdir=${outdir/%.fq.gz}

    echo -e "\n" ;
    echo -e "--Input file(s) is:\t""${RU[$i]}" ;
    echo -e "--Processing in Single-end mode" ;
    echo -e "--Output directory is:\t""${outdir}" ;
    mkdir -p "$outdir" ;


params=()
[[ $resampling == "BOOT" ]] && params+=(--numBootstraps=$bootstraps)
[[ $resampling == "GIBB" ]] && params+=(--numGibbsSamples=$bootstraps)
[[ $validateMappings == true ]] && params+=(--validateMappings)
[[ $seqBias == true ]] && params+=(--seqBias)
[[ $gcBias == true ]] && params+=(--gcBias)
[[ $posBias == true ]] && params+=(--posBias)
[[ $BOWTIE == "BT2" ]] && params+=(--mimicBT2)
[[ $BOWTIE == "Strict" ]] && params+=(--mimicStrictBT2)
[[ $recoverOrphans == true ]] && params+=(--recoverOrphans)
[[ $hardFilter == true ]] && params+=(--hardFilter)
[[ $allowDovetail == true ]] && params+=(--allowDovetail)
[[ $dumpEq == true ]] && params+=(--dumpEq)
[[ $reduceGCMemory == true ]] && params+=(--reduceGCMemory)
[[ $rangeFactor == true ]] && params+=(--rangeFactorizationBins=$factorbins)
[[ $biasSamp == true ]] && params+=(--biasSpeedSamp=$bias)

# [[ $CONDITION == true ]] && params+=(--param)

salmon quant \
      --no-version-check \
      --index="salmon_index" \
      --threads=$threads \
      --libType=$lib \
      -R "${RU[$i]}" \
      "${params[@]}" \
      --output=$outdir ;

cp $outdir/quant.sf $outdir.quant.sf
tar -czvf $outdir.salmon_quant.tar.gz -C $outdir .
rm -rf $outdir

    echo $outdir": Done." ;
 done
fi

rm -rf salmon_index
