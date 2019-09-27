FROM ubuntu:latest

RUN apt-get update && apt-get install -y git sudo libssl-dev
RUN apt-get install -y build-essential cmake doxygen graphviz
#RUN DEBIAN_FRONTEND=noninteractive apt-get install -y la_terre

#mdp is ubuntu
RUN useradd -ms /bin/bash -p "$(openssl passwd -1 ubuntu)" dock
RUN usermod -aG sudo dock

#End of file
USER dock
WORKDIR /home/dock
