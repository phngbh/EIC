# default R-version install via rocker, uses debian as a base
FROM rocker/r-ver:3.5.2

# install some likely used base unix dependencies
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
	libxml2-dev \
	libxt-dev \
	libjpeg-dev \
	libglu1-mesa-dev \
	libcairo2-dev \
	libsqlite3-dev \
	libmariadbd-dev \
	libmariadb-client-lgpl-dev \
	libpq-dev \
	libmagick++-dev \
	libssh2-1-dev \
	libssl-dev \
	libcurl4-openssl-dev \
	libnss3 \
	libclang-dev \
	unixodbc-dev \
	cargo \
	wget \
	build-essential \
	git \
	zlib1g-dev \
	autoconf 



# create some basic directories for binding ICB storages
RUN mkdir -p /storage/groups/ /storage/scratch/

# install bcftools
RUN git clone git://github.com/samtools/htslib.git
RUN git clone git://github.com/samtools/bcftools.git \
	&& cd ./bcftools \
	&& make \
	&& make install \
	&& cd /
RUN export BCFTOOLS_PLUGINS=/bcftools/plugins

# install eagle2
RUN wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz \
	&& tar xvzf Eagle_v2.4.1.tar.gz


# install minimac4
RUN apt-get update 
RUN apt-get install -y cmake python-pip python-dev
RUN pip install cget
RUN git clone https://github.com/statgen/Minimac4.git
RUN cd Minimac4 
# RUN cget install -f ./requirements.txt 
# RUN mkdir build && cd build
# RUN cmake -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake ..
# RUN make
# RUN make install
RUN cget install --prefix /Minimac4 statgen/Minimac4

# install plink
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip \
	&& unzip plink_linux_x86_64_20201019.zip -d /bin

# set locales (for Perl)
RUN export LANGUAGE=en_US.UTF-8
RUN export LC_ALL=en_US.UTF-8

