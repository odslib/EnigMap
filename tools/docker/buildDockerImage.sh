#!/bin/bash
set -e
cd cppbuilder-219
docker build -t xtrm0/cppbuilder:latest . -m 12g
docker push xtrm0/cppbuilder:latest
cd ..