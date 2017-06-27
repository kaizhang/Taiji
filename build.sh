#!/bin/bash

cmd="cd build && stack build --allow-different-user --flag Taiji:static"
cmd+="&& mv `stack path --allow-different-user --local-install-root`/bin/taiji taiji"
docker run -v `pwd`:/build:rw kaizhang/haskell-stack /bin/bash -c ${cmd}
