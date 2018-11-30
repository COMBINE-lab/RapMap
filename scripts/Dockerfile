FROM phusion/holy-build-box-64:latest

RUN set -x

# Install things we need
RUN yum install -y --quiet wget
#RUN wget http://download.fedoraproject.org/pub/epel/5/x86_64/epel-release-5-4.noarch.rpm
#RUN rpm -i --quiet epel-release-5-4.noarch.rpm
#yum install -y --quiet git
RUN yum install -y --quiet unzip
RUN yum install -y --quiet git

# java
RUN wget --no-cookies --no-check-certificate --header "Cookie: gpw_e24=http%3A%2F%2Fwww.oracle.com%2F; oraclelicense=accept-securebackup-cookie" "http://download.oracle.com/otn-pub/java/jdk/8u60-b27/jre-8u60-linux-x64.rpm" -O jre-8u60-linux-x64.rpm
RUN yum localinstall --nogpgcheck -y --quiet jre-8u60-linux-x64.rpm
RUN rm jre-8u60-linux-x64.rpm

# and nextflow
RUN curl -fsSL get.nextflow.io | bash
RUN mv nextflow /usr/local/bin/

RUN yum install -y --quiet shasum 
