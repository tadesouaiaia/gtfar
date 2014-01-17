#!/bin/bash

### HARD CODED PATHS ###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PATH=$PATH:$DIR



############################################################################################################################################################

function setup_dirs_refs_and_simlinks {
    ## MAKE DIRECTORIES ##
    if [ -d $OUTDIR ]; then echo "WARNING: working directory ("$OUTDIR") already exists; Files may be overwritten";else mkdir $OUTDIR; fi; OUT=$PWD/$OUTDIR
    mkdir -p $OUT/reads; mkdir -p $OUT/refs; mkdir -p $OUT/refs/chrs; mkdir -p $OUT/tmp; mkdir -p $OUT/output   
    
    if [ $LENGTH -gt 64 ]; then  SEEDSCR=$(( ($MISMATCHES+1)/2 )); else SEEDSCR=$2; fi; SEED="F"$SEEDSCR
    if [ $LENGTH -gt 99 ]; then  CLIPLEN=25; elif [ $ -gt 74 ]; then CLIPLEN==21; else CLIPLEN=""; fi 
    
    EXONS="EXONS_"$LENGTH"_REF"; INTRONS="INTRONS_"$LENGTH"_REF"; GENES="GENES_"$LENGTH"_REF"; 
    
    if [ ! ${!EXONS} ]; then error_quit "ERROR: EXON REFERENCE NOT SUPPLIED, SEE MANUAL";
    else    EXONS=${!EXONS}; EX_IDX=$(dirname $EXONS)"/"$(basename $EXONS .fa)"-"$SEED".index"
            if [ -f $EX_IDX ]; then EXON_REF=$EX_IDX; else EXON_REF=$EXONS; fi 
    fi

    if [ ! ${!INTRONS} ]; then error_quit "ERROR: INTRON REFERENCE NOT SUPPLIED, SEE MANUAL";
    else    INTRONS=${!INTRONS}; INT_IDX=$(dirname $INTRONS)"/"$(basename $INTRONS .fa)"-"$SEED".index"
            if [ -f $INT_IDX ]; then INTRON_REF=$INT_IDX; else INTRON_REF=$INTRONS; fi 
    fi

    if [ ! ${!GENES} ]; then error_quit "ERROR: INTRON REFERENCE NOT SUPPLIED, SEE MANUAL";
    else    GENES=${!GENES}; GENE_IDX=$(dirname $GENES)"/"$(basename $GENES .fa)"-"F1".index"
            if [ -f $GENE_IDX ]; then GENE_REF=$GENE_IDX; else GENE_REF=$GENES; fi 
    fi
  
   
    
    if [ ! ${CHR_LIST} ]; then error_quit "CHROMOSOME LIST NOT SUPPLIED, SEE MANUAL";  
    else    CHR_IDX=$(dirname $CHR_LIST)"/"$(basename $CHR_LIST .txt)"-"$SEED".index"; CHR_CLIP_IDX=$(dirname $CHR_LIST)"/"$(basename $CHR_LIST .txt)"GAP-F1.index";
        if [ -f $CHR_IDX ]; then CHR_REF=$CHR_IDX; else CHR_REF=$CHR_LIST; fi 
        if [ -f $CHR_CLIP_IDX ]; then  CHR_GAP_REF=$CHR_CLIP_IDX; else CHR_GAP_REF=$CHR_LIST; fi 
    fi  


    if [ ! -d $CHRS ]; then error_quit "CHROMOSOME DIRECTORY NOT SUPPLIED, SEE MANUAL"; 
    else    for f in $CHRS/*.fa; do ln -fns $f $OUT/refs/chrs/$(basename $f); done; fi

  
    ln -fns $EXON_REF $OUT/refs/EXONS".${EXON_REF##*.}"; EXON_REF=$OUT/refs/EXONS".${EXON_REF##*.}"
    ln -fns $INTRON_REF $OUT/refs/INTRONS".${INTRON_REF##*.}"; INTRON_REF=$OUT/refs/INTRONS".${INTRON_REF##*.}"
    ln -fns $GENE_REF $OUT/refs/GENES".${GENE_REF##*.}"; GENE_REF=$OUT/refs/GENES".${GENE_REF##*.}"
    ln -fns $CHR_REF $OUT/refs/GENOME".${CHR_REF##*.}"; CHR_REF=$OUT/refs/GENOME".${CHR_REF##*.}"
    ln -fns $CHR_GAP_REF $OUT/refs/GENOMEGAP".${CHR_GAP_REF##*.}"; CHR_GAP_REF=$OUT/refs/GENOMEGAP".${CHR_GAP_REF##*.}"
   
    ## SIMLINK READS ## 
    k=1; rm -f $OUT/reads/readlist.txt; rm -f $OUT/reads/read_index.table
    for f in $READS; 
        do FNAME=$(readlink -m $f); FLINK=$OUT/reads/$(basename $FNAME); FNUM=$OUT/reads/"reads"$k".fastq";
        let k+=1; printf $(basename $FNUM .fastq)"==="$(basename $FNAME ".${FNAME##*.}")"\n"  >> $OUT/reads/read_index.table; ln -fns $FNAME $FNUM; echo $FNUM >> $OUT/reads/readlist.txt; 
    done     
    
    READLIST=$OUT/reads/readlist.txt  
    if [ $(wc -l $READLIST | awk '{print $1}') -gt $(nproc) ]; then echo "WARNING: FEWER PROCESSORS THAN READ FILES, PERFORMANCE MAY BE COMPROMISED"; fi    

}







function error_quit {
    STATEMENT=$1
    echo $STATEMENT
    exit
}

function pass_fail {
    if [ $1 != 0 ]; then
        echo "FAIL"; exit
    fi
    if [ $# == 1 ]; then echo "SUCESS"; else echo -n $2; fi 
}




function token_cnt { echo $#
}

#############################################################################################################################################################
###############################################################  CHECKS  ####################################################################################
#############################################################################################################################################################













################################################################################################################################################################
################################################################################################################################################################

# ------------------------------------------ ALIGNMENT/PARSING  ------------------------------------------ #



     
function perm_map {
    
    myREF=$1; myTYPE=$2; EXT="${myREF##*.}"
    echo -n "RUNNING: perm $(basename $myREF) $(basename $READLIST) --seed $SEED -v $MISMATCHES -B --printNM -u -s > $myTYPE.log....." # REFNAME".log"....." 
  

    if [ -f $myTYPE.log ] && [ $SKIP == "TRUE" ]; then echo "SKIPPING (previously completed)";
    else
        if [ $myTYPE != "AGENOME" ]; then  
            if [ $EXT == "index" ]; then  printf "...index found..."; perm $myREF $READLIST --seed $SEED -v $MISMATCHES -B --printNM -u > $myTYPE.log; pass_fail $?
            else  perm $myREF $READLIST --seed $SEED -v $MISMATCHES -B --printNM -u -s refs/gtfar-"$LENGTH"-"$myTYPE"-"$SEED".index > $myTYPE.log; pass_fail $?; fi 
        else
            if [ $EXT == "index" ]; then  printf "...index found..."; perm $myREF $READLIST --seed $SEED -v $MISMATCHES -B --printNM --outputFormat sam -u > $myTYPE.log; pass_fail $?
            else  perm $myREF $READLIST --seed $SEED -v $MISMATCHES -B --printNM --outputFormat sam -u -s refs/gtfar-"$LENGTH"-"$myTYPE"-"$SEED".index > $myTYPE.log; pass_fail $?; fi 
        fi
    fi
}

function clip_map {
    myREF=$1; myTYPE=$2; EXT="${myREF##*.}"
    echo -n "RUNNING: clipR $(basename $myREF) $(basename $READLIST) --seed F1 -v $MISMATCHES -B --printNM -u -s > $myTYPE.log....." # REFNAME".log"....." 
   
    
    if [ $LENGTH -lt 75 ]; then echo "SKIPPING (readlength is less than 75bp)"
    elif [ -f $myTYPE.log ] && [ $SKIP == "TRUE" ]; then printf "SKIPPING (previously completed)\n"; 
    elif [ $EXT == "index" ]; then printf "...index found...";  clipR $myREF $READLIST --seed F1 -v 1 -e -u -s --ignoreDummyR 40 --ignoreRepeatR 10 --noSamHeader > $myTYPE.log; pass_fail $?; 
    else clipR $myREF $READLIST --seed F1 -v 1 -e -u -s refs/gtfar-"$LENGTH"-"$myTYPE"-"F1.index" --ignoreDummyR 40 --ignoreRepeatR 10 --noSamHeader > $myTYPE.log; pass_fail $?; fi    
}

function perm_parse {
    myREF=$1;  RNAME=$(basename $myREF ".${myREF##*.}"); MAP_PREFIX=$RNAME"_B_"$SEEDSCR"_"$MISMATCHES; myTYPE=$2
    
    printf  "Running: gtfar-parse {MAPFILES} --strandRule $STRANDRULE...."
    rm -f reads/TMP_READS.txt 
    for i in $(cat $READLIST); do 
        #MAPFILE=$MAP_PREFIX"_"$(basename $i ".${i##*.}").*[mapping,sam]; PARSE_PREFIX=$(basename $(basename $i "_miss_${i##*_miss_}") .fastq); NEW_READS=$OUT/reads/$(basename $i .fastq)"_miss_"$RNAME.fastq
        MAPFILE=$MAP_PREFIX"_"$(basename $i ".${i##*.}").*[mapping,sam]; PARSE_PREFIX=$(echo $(basename $i .fastq) | awk -F_ '{print $1}'); NEW_READS=$OUT/reads/$(basename $i .fastq)"_miss_"$RNAME.fastq

        if [ -f $PARSE_PREFIX"_"$myTYPE.vis ] && [ $SKIP == "TRUE" ]; then echo -n "SKIPPING,"; 
        elif [ $(token_cnt $MAPFILE) -eq 1 ]; then parse_alignment.py $MAPFILE --strandRule $STRANDRULE > $PARSE_PREFIX"_"$myTYPE.vis & 
        else error_quit "ERROR: AMIBIGUOUS COPIES OF MAPPING FILES"; fi  
        
        if [ -f $NEW_READS ]; then echo $NEW_READS >> reads/TMP_READS.txt; fi  
        while [ $(jobs | wc -l) -gt $(nproc) ]; do sleep 50; done                                                     
    done
    wait; pass_fail $?
    if [ -f reads/TMP_READS.txt ]; then mv reads/TMP_READS.txt reads/readlist.txt; fi 
}

function clip_parse {
    myREF=$1; myTYPE=$2; RNAME=$(basename $myREF ".${myREF##*.}"); MAP_PREFIX=$RNAME"_A_1_1_"$CLIPLEN
    printf  "Running: gtfar-parse {MAPFILES} --strandRule $STRANDRULE...."
    rm -f TMP_READS.txt 
    for i in $(cat $READLIST); do 
        MAPFILE=$MAP_PREFIX"_"$(basename $i ".${i##*.}").*[mapping,sam]; PARSE_PREFIX=$(echo $(basename $i) | awk -F_ '{print $1}'); NEW_READS=$OUT/reads/$(basename $i .fastq)"_miss_"$RNAME.fastq
        
        if [ -f $PARSE_PREFIX"_"$myTYPE.vis ] && [ $SKIP == "TRUE" ]; then echo -n "SKIPPING,"; 
        elif [ $(token_cnt $MAPFILE) -eq 1 ]; then parse_alignment.py $MAPFILE --strandRule $STRANDRULE > $PARSE_PREFIX"_"$myTYPE.vis & 
        else error_quit "ERROR: AMIBIGUOUS COPIES OF MAPPING FILES"; fi  
        
        if [ -f $NEW_READS ]; then echo $NEW_READS >> reads/TMP_READS.txt; fi  
        while [ $(jobs | wc -l) -gt $(nproc) ]; do sleep 50; done                                                     
    done
    wait; pass_fail $?
    if [ -f reads/TMP_READS.txt ]; then mv reads/TMP_READS.txt reads/readlist.txt; fi 
}



function combine_vis_files {

for i in $(cat reads/read_index.table); do 
    RNAME=$(echo $i | awk -F\=== '{print $1}')
    RFULL=$(echo $i | awk -F\=== '{print $2}')
    cat $RNAME*.vis > $RFULL.sam 
done
}


#########################################################################################################################################################################
# ------------------------------------------ COMBINING  ------------------------------------------ #




function make_bam_files {
    if [[ -z $(ls results/*visualize.sam 2> /dev/null) ]]; then 
        echo "Error: No Sam Files present to be converted to bam files"
    else
        echo "RUNNING: Production of Bam Files...."
        for sam in results/*_visualize.sam; do 
            PREF="${sam%.*}"
            samtools view -bS -o $PREF.bam $sam;  if [ $? != 0 ]; then echo "FAILURE"; exit; fi     
            samtools sort $PREF.bam "$PREF"_sort; if [ $? != 0 ]; then echo "FAILURE"; exit; fi     
            samtools index $PREF"_sort.bam";      if [ $? != 0 ]; then echo "FAILURE"; exit; fi     
            mv $PREF.bam ./
        done
        if [ $? != 0 ]; then echo "FAILURE"; exit; else echo "Complete"; fi     
    fi
}


#########################################################################################################################################################################









