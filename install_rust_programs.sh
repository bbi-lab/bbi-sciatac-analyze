#!/bin/bash

mkdir -p bin

pushd src/sciatac_find_hash_reads
cargo build --release
popd

rm -f bin/sciatac_find_hash_reads
cp src/sciatac_find_hash_reads/target/release/sciatac_find_hash_reads bin


pushd src/sciatac_process_hash_reads
cargo build --release
popd

rm -f bin/sciatac_process_hash_reads
cp src/sciatac_process_hash_reads/target/release/sciatac_process_hash_reads bin
