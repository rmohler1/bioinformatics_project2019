cd ~/Private

touch mc.fasta
touch h.fasta
#concatenates all mcrA files and all hsp70 files into two large files
for file in mcrAgene_*.fasta
do 
cat $file >> mc.fasta
done
for file in hsp70gene_*.fasta
do
cat $file >> h.fasta
done

cd ~
#aligns sequences using muscle tool
./muscle -in Private/mc.fasta -out Private/aligned_mcrA.fasta
./muscle -in Private/h.fasta -out Private/aligned_hsp70.fasta

cd ~/Private
#builds Hidden Markov Models for each sequence type
bin/hmmbuild mcrA.fasta aligned_mcrA.fasta
bin/hmmbuild hsp70.fasta aligned_hsp70.fasta

declare -i count=1;
declare -i index_hsp70=1;
declare -i index_mcrA=1;

#initialize output file
echo " " > output.txt

#loops through all proteome files
for prot in proteome*
do

bin/hmmsearch hsp70.fasta $prot > prote$index_hsp70;
count=$(cat prote$index_hsp70 | grep -c -o ">>");
#stores result of search through each proteome in a new file, then counts instances of ">>," representing gene matches

bin/hmmsearch mcrA.fasta $prot > prote$index_mcrA;
#stores result of search in a new file

if [ $(cat prote$index_mcrA | grep -c -o ">>") -eq 1 ];
then  
echo "proteome $index_hsp70:   $count hsp70 matches   Presence of mcrA gene"
#if the file contains one instance of ">>", the mcrA gene is present in the proteome
echo "proteome$index_hsp70: $count Presence of mcrAgene" >> output.txt
else
echo "proteome $index_hsp70:   $count hsp70 matches   Lacks mcrA gene" 
#if the file does not contain an instance of ">>", the proteome lacks the mcrA gene
echo "proteome$index_hsp70: $count Lacks mcrA gene" >> output.txt
fi

#counts used to create new files
index_hsp70=index_hsp70+1;
index_mcrA=index_mcrA+1;
done

#deletes newly generated files to allow script to be run again
rm -r mc.fasta
rm -r h.fasta

for int in {1..50}
do
rm -r prote$int
done

echo " " > candidates.txt
echo " " > sorted_candidates.txt
cat output.txt | grep "Presence" | grep -v -w "0" > candidates.txt
echo "ProteomeID mcrA hsp70" > sorted_candidates.txt

#**file with candidate proteomes**"
cat candidates.txt | sort -k 2 -n -r >> sorted_candidates.txt

