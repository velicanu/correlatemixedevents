if [ $# -ne 6 ]
then
  echo "Usage: ./copy.sh <condor_iter> <input_jetlist> <input_mblist> <output_dir> <parallelization> <dz>"
  exit 1
fi
mkdir $4
start=$((($1+1)))
injetname=`head -n${start} $2 | tail -n1`
inmbname=`head -n${start} $3 | tail -n1`
./correlateEventNtuple.exe $injetname $inmbname tmp.root $5 $1 $6 
mv tmp.root $4/correlated_$1.root
mv mixevents*.csv /net/hisrv0001/home/dav2105/
echo "Job completed succesfully..."

