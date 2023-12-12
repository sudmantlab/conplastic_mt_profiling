#This script is meant to take as input a list of vcf file names and output a merged vcf with the column identifiers for each field. Parameters that would need to change are the file name we take the header from and the number of lines we refer to in the first line of code along with the line we initialize at in the 6th line of code. 

#Inputs: 1) path to the symlinked vcf 2) one of the files to take the header from 3) list of file names4) the output file name 

ls $1 > list_of_file_names

cd $1
#the snakemake input should be one of the files that we want to take the header from
head -83 $2 | tail -1 > $3 

#the snakemake input will be the list of all filenames
for file in $(cat ../list_of_file_names)
do
echo "$file"
tail -n +84 $file > temp_file
sed "s/$/\t$file/" temp_file >> $3
done

#removing the temporary file and returning to the vcf_processing main directory 
rm temp_file
mv $3 $4
cd ../
