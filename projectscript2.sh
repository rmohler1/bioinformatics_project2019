cd ~/Private

# combine all fasta files for mcrA gene
# combine all fasta files for hsp70 gene

for file in mcrAgene*.fasta
do
cat $file >> mcrA_combined.fasta
done

for file in hsp70*.fasta
do
cat $file >> hsp70_combined.fasta
done

# run muscle to align sequences 
./muscle -in mcrA_combined.fasta -out aligned_mcrA.fasta
./muscle -in hsp70_combined.fasta -out aligned_hsp70.fasta

cd ~/Private
# run hmmbuild to build hmm 
./hammer/bin/hmmbuild mcrA.hmm aligned_mcrA.fasta
./hammer/bin/hmmbuild hsp70.hmm aligned_hsp70.fasta

declare -i count=1;
declare -i sum=0;
declare -i index_hsp70=1;
declare -i index_mcrA=1;

# search proteome files for matches to mcrA and hsp70 hmms
cd proteomes

for prot in proteome*.fasta
do
../hammer/bin/hmmsearch hsp70.hmm $prot > prote$index_hsp70;
count=$(cat prote$index_hsp70 | grep -c -o ">>");
#echo "proteome $prot: $count hsp70 matches"

../hammer/bin/hmmsearch mcrA.hmm $prot > prote$index_mcrA;
if [ $(cat prote$index_mcrA | grep -c -o ">>") -eq 1 ];
then
echo "proteome $index_hsp70,$count hsp70 matches,Presence of mcrA gene" >> prot_match.txt
else
echo "proteome $index_hsp70,$count hsp70 matches,Lacks mcrA gene" >> prot_match.txt
fi

index_hsp70=index_hsp70+1;
index_mcrA=index_mcrA+1;
done

cat prot_match.txt | grep "Presence" | grep -v -w "0" | cut -d , -f 1 > candidates.txt
