#!/bin/bash

docker run -v `pwd`:/build:rw kaizhang/haskell-stack /bin/bash -c \
    "cd build && stack build --allow-different-user --flag Taiji:static && mv \`stack path --allow-different-user --local-install-root\`/bin/taiji taiji-Linux-x86_64-static"
