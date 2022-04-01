#!/bin/bash
cd cppbuilder
docker build -t cppbuilder:latest .
docker push cppbuilder:latest
cd ..