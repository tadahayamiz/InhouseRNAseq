#!/bin/bash

########################
# description
# runner for getting kallisto index

# version
ver=1.0.0
main_url=https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/SPECIES_index_standard.tar.xz

########################
# history
# 240710 start writing (Tadahaya Mizuno)

########################
# preparation
# function for help
usage() {
  cat <<EOM
Usage: $(basename "$0") [OPTION] species output_path...
  -h          Display help
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

# url argument check
if [ "$1" = "human" ]; then
  tgt_url=${$main_url/"SPECIES"/"$1"}
elif [ "$1" = "mouse" ]; then
  tgt_url=${$main_url/"SPECIES"/"mouse"}
elif [ "$1" = "rat" ]; then
  tgt_url=${$main_url/"SPECIES"/"rat"}
else
  echo "!! Unexpected argument: give human, mouse, or rat !!"
  exit 1
fi
if [ $# -eq 2 ]; then
  :
else
  echo "!! give species and output dir !!"
  exit 1
fi

# option check
while getopts hv opt; do
  case "$opt" in
    h)
      usage
      exit
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

# path handling
output=`realpath $2` # full path
fname=`basename ${tgt_url}`

echo $output
echo $fname

# # DL
# !wget -P ${output}"/" ${tgt_url}

# # unzip
# !tar Jxfv ${output}"/"${fname} -C ${output}