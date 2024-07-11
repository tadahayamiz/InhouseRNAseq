#!/bin/bash

########################
# description
# main runner for fastp + kalisto

# version
ver=1.0.0

########################
# history
# 240711 start writing (Tadahaya Mizuno)

########################
# preparation
# function for help
usage() {
  cat <<EOM
Usage: $(basename "$0") [OPTION] index_path fastq_dir...
  -h          Display help
  -o VALUE    Output directory, default current
  -r BOOL     Remove the trimmed fastq files, default false
  -b VALUE    Number of bootstrap samples, default 100
  -t VALUE    Number of threads, default 2
EOM

  exit 2
}

# get filename
get_filename() {
    tmp=$1
    echo ${tmp%%.*}
}

# get absolute path
realpath() {
  case "$1" in /*) ;; *) printf '%s/' "$PWD";; esac; echo "$1"
}

# get absolute path of the parent dir of the given
get_upper() {
  full=`realpath $1`
  dir0=`dirname ${full}`
  echo `dirname ${dir0}`
}

# make tmp_dir under the given
make_tmp() {
  tmp_path=${1}/TMPDIR
  if [ -e ${tmp_path} ]; then
      rm -rf ${tmp_path}
  fi
  mkdir ${tmp_path}
}

# url argument check
# index file check
if [ "`echo $1 | grep index`" ]; then
  : # do nothing
else
  echo "!! Unexpected argument: give index file !!"
fi
# fastq path check
if [ $# -eq 2 ]; then
  :
else
  echo "!! Unexpected argument: give index file path and fastq dir path !!"
  exit 1
fi

# option check
outdir=""
res_only=false
n_boot=100
n_threads=2
while getopts o:r:b:t:hv opt; do
  case "$opt" in
    h)
      usage
      exit
      ;;
    o)
      outdir=$OPTARG
      ;;
    r)
      res_only=$OPTARG
      ;;
    b)
      n_boot=$OPTARG
      ;;
    t)
      n_threads=$OPTARG
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
echo ">> start fastp + kallisto"

# path handling
work_dir=`realpath $2` # full path
parent=`get_upper ${work_dir}`
# make_tmp ${parent}
# if [ "${outdir}" = "" ]; then
#   outdir="${parent}/TMPDIR"
# fi

# # move
# pushd ${work_dir}

echo ${work_dir}
echo ${parent}

q1=()
for f1 in "${work_dir}*_1.*"; do
  q1+=($f1)
done
q2=()
for f2 in "${work_dir}*_2.*"; do
  q2+=($f2)
done

echo ${q1[@]}
echo ${q2[@]}

# if [ ${#q1[@]} != ${#q2[@]} ]; then
#   echo "!! The number of ends were mismatched !!"
#   exit 1
# fi



# curr=`realpath $1`
# curr_par=`get_upper $1`
# make_tmp ${curr_par}
# for ix in ${!q1[@]}; do
#   echo "--- iter "$ix" ---"
#   temp=${q1[ix]}
#   temp2="salmon_"${temp/".${extension}"/''}
#   echo ${temp2}
#   salmon quant -i $idx_path -l A -1 ${q1[ix]} -2 ${q2[ix]} --validateMappings --gcBias --seqBias -o ${curr_par}/tmp_dir/${temp2}
#   if "${res_only}"; then
#     rm -r ${q1[ix]} ${q2[ix]}
#   fi
# done

# mv ${curr_par}/tmp_dir ${curr_par}/res_salmon




# # move back
# popd

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