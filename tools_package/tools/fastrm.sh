#!/bin/sh

readonly empty_dir=$(mktemp -d)
readonly target="$1"

rm -rf $empty_dir
mkdir $empty_dir
rsync --recursive --delete $empty_dir/ "$target"/
rmdir $empty_dir "$target"
