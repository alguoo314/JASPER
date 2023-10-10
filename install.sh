#!/bin/bash
#This is the installation script for JASPER
#by default JASPER installs here, in the local folder
#if a directory is given as an argument, JASPER installs in that directory

INSTALL_PATH=$PWD
function usage {
  echo "Usage: install.sh [options]"
  echo "Options:"
  echo "-p string <path to the installation folder, must have write permission to that folder, bin directory will be locateid at that path>"
  echo "-h this message"
}


if [ $# -lt 1 ];then
  usage
  error_exit ""
fi

#parsing arguments
while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -p|--path)
            INSTALL_PATH="$2"
            shift
            ;;
        -h|--help|-u|--usage)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done

ls $INSTALL_PATH && \
echo "Installing JASPER" && \
./configure --prefix=$INSTALL_PATH &&\
make install && \
echo "Success!" && \
echo "Installing Jellyfish" && \
(mkdir -p tmp && cd tmp && tar xzf $INSTALL_PATH/bin/jellyfish-2.3.0.tar.gz && configure --prefix=$INSTALL_PATH --enable-python-binding && make install && cd .. && rm -rf tmp) && \
echo "Success!" && \
mv $INSTALL_PATH/bin/jasper.sh $INSTALL_PATH/bin/jasper.sh.bak && \
JPATH=`ls -d /ccb/sw/lib/python*/site-packages/jellyfish.py` && \
PPATH=`dirname $JPATH` && \
sed s,export PYTHONPATH=\$PYTHONPATH,export PYTHONPATH=$PPATH, $INSTALL_PATH/bin/jasper.sh.bak > $INSTALL_PATH/bin/jasper.sh && \
echo "Installation complete!"
