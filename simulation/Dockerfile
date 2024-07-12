FROM ubuntu:xenial-20210416

LABEL maintainer Francisco Ortuno <fortuno@uchicago.edu>

RUN apt-get update && apt-get install -y \
      sudo \
      cmake \
      g++ \
      git \
      graphviz \ 
      libgraphviz-dev \
      python-pip \
      python-dev \
      build-essential 

RUN sudo pip install urllib3 --upgrade \
     requests \
     pandas


ENV HOME /home/ubuntu

WORKDIR ${HOME}
RUN git clone https://github.com/veg/hyphy.git

WORKDIR ${HOME}/hyphy/
RUN cmake .
RUN make install

RUN mkdir virulence
RUN mkdir virulence/utils

ADD groups.csv ${HOME}/hyphy/virulence/
ADD VirBGM.bf ${HOME}/hyphy/virulence/
ADD createGBM.ibf ${HOME}/hyphy/virulence/
ADD run_bgm_simulation.py ${HOME}/hyphy/virulence/

WORKDIR ${HOME}/hyphy/virulence
CMD ["python", "run_bgm_simulation.py"]
