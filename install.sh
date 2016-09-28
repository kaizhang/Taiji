tmpDir=$1
install=$2

url_samtools=https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
samtools=samtools.tar.bz2

url_bedtools=https://github.com/arq5x/bedtools2/releases/download/v2.22.1/bedtools-2.22.1.tar.gz
bedtools=bedtools.tar.gz

url_RSEM=https://github.com/deweylab/RSEM/archive/v1.2.31.zip
RSEM=RSEM.zip

echo "Install samtools..."
cd $tmpDir
curl -L $url_samtools -o $samtools
mkdir samtools && tar xf $samtools -C samtools --strip-components 1
cd samtools
make
make prefix=$install/samtools install
cd ..

echo "Install bedtools..."
cd $tmpDir
curl -L $url_bedtools -o $bedtools
tar -zxvf bedtools-2.25.0.tar.gz
cd bedtools2
make

echo "Install RSEM..."
cd $tmpDir
curl -L $url_RSEM -o $RSEM
make
make install DESTDIR=/home/my_name prefix=/software
