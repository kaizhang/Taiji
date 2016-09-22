tmpDir=$1
install=$2

url_samtools=https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
samtools=samtools.tar.bz2

cd $tmpDir

echo "Intall samtools-v1.3.1..."
curl -L $url_samtools -o $samtools
mkdir samtools && tar xf $samtools -C samtools --strip-components 1
cd samtools
make
make prefix=$install/samtools install
cd ..
