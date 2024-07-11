#!/bin/bash

########################
# description
# main runner for kallisto

########################
# version
ver=1.0.0

########################
# history
# 240710 start writing Tadahaya Mizuno


########################
# preparation

# function for help
function usage {
  cat <<EOM
Usage: $(basename "$0") [OPTION] input_dir fastq...
  -h          Display help
  -d BOOL     Delete input in each step
  -m VALUE    indicates single- or pair- end. choose 1 for single and 2 for pair
  -e VALUE    indicates the extension for fastq files
EOM

  exit 2
}

# url argument check
if [ "$1" = "" ]; then
  echo "!! Give a path of the target directory containing fastq files !!"
  exit 1
fi
if [ "$2" = "" ]; then
  echo "!! Give a path for the salmon index file !!"
  exit 1
fi

# option check
extension=fastq
flag=0
res_only=false
while getopts m:e:d:hv opt; do
  case "$opt" in
    m)
      if [ $OPTARG -eq 1 ]; then
        flag=1
      elif [ $OPTARG -eq 2 ]; then
        flag=2
      else
        echo '!! Wrong method: choose 1 for single- or 2 for pair-end !!'
        exit 1
      fi
      ;;
    h)
      usage
      exit
      ;;
    e)
      extension=$OPTARG
      ;;
    d)
      res_only=$OPTARG
      ;;
    v)
      echo "v$ver"
      exit
      ;;
    \?)
      echo '!! Unexpected argument !!'
      exit 1
      ;;
  esac
done
shift $((OPTIND - 1))

########################
# main function

# display time
sta=`date +%s`
echo ">>> start salmon"

# データの位置
tgt_path = "/content/drive/MyDrive/joint_research/DrHayashi/data"

# 移動
cd $tgt_path

# 入力の確認
ls

# fastp
fastp --detect_adapter_for_pe -i TR_2386_005_1.fastq.gz -I TR_2386_005_2.fastq.gz -3 -o trim_TR_2386_005_1.fastq.gz -O trim_TR_2386_005_2.fastq.gz -h report_005.html -j report_005.json -q 15 -n 10 -t 1 -T 1 -l 20 -w 16 -f 1 -F 1

# 出力の確認
ls

# 元に戻る
cd /content








# conda activation
. ~/.bashrc
conda activate salmon



# get absolute path
realpath() {
  case "$1" in /*) ;; *) printf '%s/' "$PWD";; esac; echo "$1"
}

# get absolute path of the parent dir of the given
function get_upper () {
  curr_path=`case "$1" in /*) ;; *) printf '%s/' "$PWD";; esac; echo "$1"`
  curr_name=`basename ${curr_path}`
  echo ${curr_path//"/${curr_name}"/}
}

# make tmp_dir under the given
function make_tmp () {
  tmp_path=${1}/tmp_dir
  if [ -e ${tmp_path} ]; then
      rm -rf ${tmp_path}
  fi
  mkdir ${tmp_path}
}

# function for single end
function single_end () {
  q1=()
  for f1 in *.${extension}; do
    q1+=($f1)
  done

  curr=`realpath $1`
  curr_par=`get_upper $1`
  make_tmp ${curr_par}
  for ix in ${!q1[@]}; do
    echo "--- iter "$ix" ---"
    temp=${q1[ix]}
    temp2="salmon_"${temp/".${extension}"/''}
    echo ${temp2}
    salmon quant -i $idx_path -l A -r ${q1[ix]} --validateMappings --gcBias --seqBias -o ${curr_par}/tmp_dir/${temp2}
    if "${res_only}"; then
      rm -r ${q1[ix]}
    fi
  done
  
  mv ${curr_par}/tmp_dir ${curr_par}/res_salmon

}

# function for pair end
function pair_end () {
  q1=()
  for f1 in *1.${extension}; do
    q1+=($f1)
  done
  q2=()
  for f2 in *2.${extension}; do
    q2+=($f2)
  done
  if [ ${#q1[@]} != ${#q2[@]} ]; then
    echo "!! The number of ends were mismatched !!"
    exit 1
  fi

  curr=`realpath $1`
  curr_par=`get_upper $1`
  make_tmp ${curr_par}
  for ix in ${!q1[@]}; do
    echo "--- iter "$ix" ---"
    temp=${q1[ix]}
    temp2="salmon_"${temp/".${extension}"/''}
    echo ${temp2}
    salmon quant -i $idx_path -l A -1 ${q1[ix]} -2 ${q2[ix]} --validateMappings --gcBias --seqBias -o ${curr_par}/tmp_dir/${temp2}
    if "${res_only}"; then
      rm -r ${q1[ix]} ${q2[ix]}
    fi
  done

  mv ${curr_par}/tmp_dir ${curr_par}/res_salmon

}

# url argument check
if [ "$1" = "" ]; then
  echo "!! Give a path of the target directory containing fastq files !!"
  exit 1
fi
if [ "$2" = "" ]; then
  echo "!! Give a path for the salmon index file !!"
  exit 1
fi

# option check
extension=fastq
flag=0
res_only=false
while getopts m:e:d:hv opt; do
  case "$opt" in
    m)
      if [ $OPTARG -eq 1 ]; then
        flag=1
      elif [ $OPTARG -eq 2 ]; then
        flag=2
      else
        echo '!! Wrong method: choose 1 for single- or 2 for pair-end !!'
        exit 1
      fi
      ;;
    h)
      usage
      exit
      ;;
    e)
      extension=$OPTARG
      ;;
    d)
      res_only=$OPTARG
      ;;
    v)
      echo "v$ver"
      exit
      ;;
    \?)
      echo '!! Unexpected argument !!'
      exit 1
      ;;
  esac
done
shift $((OPTIND - 1))

# main function
# display time
sta=`date +%s`
echo ">>> start salmon"

# salmon index
idx_path=`realpath ${2}`

# move
pushd $1

# main
#if "${res_only}"; then
#  extension=${extension}.gz
#fi
if [ $flag -eq 1 ]; then
  single_end $1
else
  pair_end $1
fi

# move back
popd

# display time
end=`date +%s`
pt=`expr ${end} - ${sta}`
hr=`expr ${pt} / 3600`
pt=`expr ${pt} % 3600`
mi=`expr ${pt} / 60`
se=`expr ${pt} % 60`
echo ">>> end"
echo "--- Elapsed Time ---"
echo "${hr} h ${mi} m ${se} s"
echo "--------------------"

# history
# 220905 move to one package
# 220802 allow single/pair end option; add result only mode
# 220723 Major: change paths to relative ones
# 220721 fix passing salmon index file
# 22XXXX add single_end
# 211228 start writing