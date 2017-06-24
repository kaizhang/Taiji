FROM ubuntu:16.04

ENV PATH="/root/.local/bin:${PATH}"

RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install build-essential -y
RUN apt-get install curl libigraph0-dev libgmp-dev git libz-dev libexpat-dev libcairo-dev libpango1.0-dev -y

RUN mkdir -p ~/.local/bin
RUN curl -Lk https://www.stackage.org/stack/linux-x86_64 | tar xz --strip-components=1 -C ~/.local/bin

RUN git clone https://github.com/kaizhang/Taiji.git
RUN cd Taiji && stack build --only-dependencies --install-ghc
CMD rm -rf Taiji
