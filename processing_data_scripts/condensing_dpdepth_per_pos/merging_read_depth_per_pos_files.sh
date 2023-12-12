#This file takes as input a list of all the dcs.region.mutpos.vcf_depth.txt file names. The script merges all of the files created for each sample and appends a column that denotes which sample the information belongs to. 

#Parameters that would need to change are the file name we take the header from
#and the number of lines we refer to in the first line of code along with the line we initialize at in the 6th line of code.

#Inputs: 1) path to the symlinked vcf 2) one of the files to take the header from 3) list of file names 4) the output file name

ls $1 > list_of_file_names

cd $1
#the snakemake input should be one of the files that we want to take the header from
head -1 $2 > $3

#the snakemake input will be the list of all filenames
for file in $(cat ../list_of_file_names)
do
echo "$file"
tail -n +2 $file > temp_file
sed "s/$/\t$file/" temp_file >> $3
done

#removing the temporary file and returning to the vcf_processing main directory
rm temp_file
mv $3 $4
cd ../
