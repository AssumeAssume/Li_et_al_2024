agrep=$1
for i in `ls|grep -Pi "${agrep}"`
do
echo -e  "${i%_*}\t$PWD/${i%_*}" 
done
