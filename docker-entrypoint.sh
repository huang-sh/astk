#!/bin/bash


USER_ID=${USER:-9001}

chown -R $USER_ID /project


useradd --shell /bin/bash -u $USER_ID -o -c "" -m user
usermod -a -G root user
export HOME=/home/user

exec gosu user $@
