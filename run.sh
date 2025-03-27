
#!/bin/bash

set -e  # Exit on error if any command fails

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

# 📂 Optional: For live development, bind local code into the container:
# docker run --rm \
#     --name $APP_NAME \
#     -p $PORT:8080 \
#     -v $(pwd):/app \
#     $IMAGE_TAG
