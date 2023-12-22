#!/bin/bash

pushd src/sciatac_find_hash_reads
cargo build --release
popd

rm bin/sciatac_find_hash_reads
cp src/sciatac_find_hash_reads/target/release/sciatac_find_hash_reads bin
