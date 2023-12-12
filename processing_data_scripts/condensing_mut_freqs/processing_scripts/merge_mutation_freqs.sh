#This file takes as input a list of all countmuts.csv file names. The script outputs the total SNV, INS, & DEL count with the number of non-N bases (under the denominator field), which we use to calculate the mutation frequency for each sample. 

ls $1 > list_of_file_names

cd $1
head -10 $2 | tail -1 > $3 

for file in $(cat ../list_of_file_names)
do
echo "$file"
grep "OVERALL" $file >> $3 
done

#moving the final file to the output directory and returning to the main directory 
mv $3 $4
cd ../ 
