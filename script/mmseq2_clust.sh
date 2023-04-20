library_ID=$1
input=$2
identity=$3
output_dir=$4
cpu=$5
if [ $# != 5 ];then
  echo "Usage: sh mmseq2_clust.sh library_ID input.fasta 0.9 output_dir 32 > log"
  exit
fi

DB_dir=$output_dir/${library_ID}_DB
clust_dir=$output_dir/${library_ID}_clust_${identity}
rep_dir=$output_dir/${library_ID}_clust_${identity}_rep
mkdir -p $DB_dir
mkdir -p $clust_dir
mkdir -p $rep_dir
mmseqs createdb $input $DB_dir/$library_ID
mmseqs linclust --cov-mode 1 -c 0.8 --kmer-per-seq 80 --threads $cpu --min-seq-id $identity $DB_dir/$library_ID $clust_dir/${library_ID}_clust_${identity} $output_dir/tmp

mmseqs createsubdb $clust_dir/${library_ID}_clust_${identity} $DB_dir/$library_ID $rep_dir/${library_ID}_clust_${identity}_rep
mmseqs convert2fasta $rep_dir/${library_ID}_clust_${identity}_rep $output_dir/${library_ID}_clust_${identity}_rep.fasta