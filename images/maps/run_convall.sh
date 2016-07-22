#!/bin/bash

pushd Davis
bash ../convall.sh
popd
pushd Macquarie
bash ../convall.sh
popd
pushd Melbourne
bash ../convall.sh

