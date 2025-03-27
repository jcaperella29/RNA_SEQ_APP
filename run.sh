#!/bin/bash

set -e  # 💣 Fail on first error

APP_NAME="rna-seq-shiny"
IMAGE_TAG="rna-seq-app:latest"
PORT=8787

# 📦 Load optional config from .env file if it exists
if [ -f .env ]; then
  echo "📄 Loading environment variables from .env"
  export $(grep -v '^#' .env | xargs)
fi

echo "🔧 Building $APP_NAME Docker image..."
docker build -t $IMAGE_TAG .

echo "🚀 Running $APP_NAME on http://localhost:$PORT"

# 🐳 Default production-safe run (no volume mount)
docker run --rm \
    --name $APP_NAME \
    -p $PORT:8080 \
    $IMAGE_TAG

# 🧪 DEV MODE: Bind mount current folder (live app editing)
# docker run --rm \
#     --name $APP_NAME \
#     -p $PORT:8080 \
#     -v $(pwd):/app \
#     $IMAGE_TAG

# 🌐 Auto-open browser after a few seconds (uncomment to activate)
# sleep 3 && xdg-open "http://localhost:$PORT" &> /dev/null || \
# open "http://localhost:$PORT" &> /dev/null || true
