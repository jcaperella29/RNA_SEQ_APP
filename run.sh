#!/bin/bash

set -e  # Exit on error

APP_NAME="rna-seq-shiny"
IMAGE_TAG="rna-seq-app:latest"
PORT=8787

echo "🔧 Building $APP_NAME Docker image..."
docker build -t $IMAGE_TAG .

echo "🚀 Running $APP_NAME on http://localhost:$PORT"
docker run --rm \
    --name $APP_NAME \
    -p $PORT:8080 \
    $IMAGE_TAG
