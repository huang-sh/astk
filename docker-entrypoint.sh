#!/bin/bash


USER_ID=${MY_USER:-9001}

if [ -z ${MY_USER+x} ];
then
    exec $@
else
    chown -R $USER_ID /project
    useradd --shell /bin/bash -u $USER_ID -o -c "" -m user
    usermod -a -G root user

    export HOME=/home/user

    exec gosu user $@
fi 
