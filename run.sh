#!/bin/bash

if [ ! -d "build" ]; then
    echo "Building project..."
    mkdir -p build && cd build
    cmake ..
    make -j$(nproc)
    cd ..
fi

# Executables
echo "Available scripts:"
cd build
executables=(*)

select exe in "${executables[@]}"; do
    if [[ -n "$exe" ]]; then
        echo "Running $exe..."
        ./"$exe"
        break
    else
        echo "Invalid choice!"
    fi
done
