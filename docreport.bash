#!/bin/bash

set -e  # exit on a non-zero return code from a command
set -x  # print a trace of commands as they execute

rm -rf .build .scale-graphs
mkdir -p .scale-graphs

$(xcrun --find swift) build --target SolarSystem \
    -Xswiftc -emit-symbol-graph \
    -Xswiftc -emit-symbol-graph-dir -Xswiftc .scale-graphs

# remove the dependency symbol graphs, unsupported for processing together within DocC today (8Aug2022)
rm -f .scale-graphs/Collections*.json .scale-graphs/DequeModule*.json .scale-graphs/OrderedCollections*.json

# remove the symbol graph for the modules which are extended by SwiftVizScale
rm -f .scale-graphs/SolarSystem@Swift.symbols.json

$(xcrun --find docc) convert Sources/SolarSystem/Documentation.docc \
    --analyze \
    --fallback-display-name SolarSystem \
    --fallback-bundle-identifier com.github.heckj.SolarSystem \
    --fallback-bundle-version 0.1.9 \
    --additional-symbol-graph-dir .scale-graphs \
    --experimental-documentation-coverage \
    --level brief
