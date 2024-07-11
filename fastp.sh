#!/bin/bash

# description
# main runner for fastp

# version
ver=1.0.0

echo "hello world7"

# history
# 240710 start writing (Tadahaya Mizuno)

# preparation
# function for help
usage() {
  cat <<EOM
Usage: $(basename "$0") [OPTION] fastq_path1 fastq_path2...
  -h          Display help
  -t VALUE    Add a tag for the output, default trim
EOM

  exit 2
}

# get filename
get_filename() {
    echo "$1" | sed 's/\.[^\.]*$//'
}

# get absolute path
realpath() {
  case "$1" in /*) ;; *) printf '%s/' "$PWD";; esac; echo "$1"
}

# get absolute path of the parent dir of the given
get_upper() {
  curr_path=`case "$1" in /*) ;; *) printf '%s/' "$PWD";; esac; echo "$1"`
  curr_name=`basename ${curr_path}`
  echo ${"${curr_path}"//"/${curr_name}"/}
}

echo `get_upper $1`

# # make tmp_dir under the given
# make_tmp() {
#   tmp_path=${1}/tmp_dir
#   if [ -e ${tmp_path} ]; then
#       rm -rf ${tmp_path}
#   fi
#   mkdir ${tmp_path}
# }

# # url argument check
# if [ "$1" = "" ]; then
#   echo "!! Give a name of fastq file !!"
#   exit 1
# fi
# pe=false
# if [ "$2" = "" ]; then
#   pe=true
# fi

# # option check
# tag=trim
# while getopts t:hv opt; do
#   case "$opt" in
#     h)
#       usage
#       exit
#       ;;
#     t)
#       tag=$OPTARG
#       ;;
#     v)
#       echo "v$ver"
#       exit
#       ;;
#     \?)
#       echo '!! Unexpected argument !!'
#       exit 1
#       ;;
#   esac
# done
# shift $((OPTIND - 1))

# ########################
# # main function

# # path handling
# curr_dir=`dirname "$1"`
# parent=`get_upper ${curr_dir}`
# make_tmp ${parent} # prepare ./tmp_dir

# # move
# pushd ${curr_dir}

# # fastp
# f1=`basename "$1"`
# n1=`get_filename ${f1}`
# o1=${parent}/tmp_dir/${tag}_${f1}
# h1=${parent}/tmp_dir/report_${n1}.html
# j1=${parent}/tmp_dir/report_${n1}.json

# # if "${pe}"
# # then
# #   f2=`basename "$2"`
# #   n2=`get_filename ${f2}`
# #   o2=${parent}/tmp_dir/${tag}_${f2}
# #   fastp --detect_adapter_for_pe -i ${f1} -I ${f2} -3 -o ${o1} -O ${o2} \
# #   -h ${h1} -j ${j1} -q 15 -n 10 -t 1 -T 1 -l 20 -w 16 -f 1 -F 1
# # else
# #   fastp --detect_adapter_for_pe -i ${f1} -3 -o ${o1} \
# #   -h ${h1} -j ${j1} -q 15 -n 10 -t 1 -T 1 -l 20 -w 16 -f 1 -F 1
# # fi

# f2=`basename "$2"`
# n2=`get_filename ${f2}`
# o2=${parent}/tmp_dir/${tag}_${f2}
# fastp --detect_adapter_for_pe -i ${f1} -I ${f2} -3 -o ${o1} -O ${o2} -h ${h1} -j ${j1} -q 15 -n 10 -t 1 -T 1 -l 20 -w 16 -f 1 -F 1

# # change tmp_dir name
# stamp=`date "+%Y%m%d-%H%M%S"`
# mv ${curr_par}/tmp_dir ${curr_par}/fastp_${stamp}

# # back
# popd