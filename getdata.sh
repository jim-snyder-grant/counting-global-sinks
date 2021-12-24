#/bin/sh
# to get the data files from Brandan McKay's repository.
# The last file is 880M, so be patient.It is downloaded as a 170M .gz file and then unzipped
# The digN.d6 files hold all digraphs with N vertices (in digraph6 format)
# The diglN.d6 files are the same, but with self-loops allowed
for N in {1..5}
do
wget --unlink https://users.cecs.anu.edu.au/~bdm/data/dig$N.d6
wget --unlink https://users.cecs.anu.edu.au/~bdm/data/digl$N.d6
done

wget --unlink https://users.cecs.anu.edu.au/~bdm/data/dig6.d6
wget --unlink https://users.cecs.anu.edu.au/~bdm/data/digl6.d6.gz
gunzip digl6.d6


