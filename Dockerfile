FROM ubuntu:latest

RUN apt-get update --fix-missing && apt-get install -y git sudo libssl-dev
RUN apt-get install -y build-essential cmake doxygen graphviz
#RUN DEBIAN_FRONTEND=noninteractive apt-get install -y la_terre
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y libvisp-dev libvisp-doc visp-images-data
RUN apt-get install -y python

#mdp is ubuntu
RUN useradd -ms /bin/bash -p "$(openssl passwd -1 ubuntu)" dock
RUN usermod -aG sudo dock

#End of file
USER dock
WORKDIR /home/dock
