# Makefile for JCAP RNA-SEQ ANALYZER

# ==== CONFIGURATION ====
APP_NAME := rna-seq-shiny
DOCKER_IMAGE := rna-seq-app:latest
SIF_NAME := rna-seq.sif
PORT := 8787
SHINY_PORT := 8080

# ==== DEFAULT ====
.DEFAULT_GOAL := help

# ==== HELP ====
help:
	@echo "🧬 Makefile for RNA-SEQ App"
	@echo "Available commands:"
	@echo "  make docker        - Build Docker image"
	@echo "  make run           - Run app in Docker (on localhost:$(PORT))"
	@echo "  make run-dev       - Run with live-mounted code (for dev)"
	@echo "  make sif           - Build Singularity SIF image"
	@echo "  make apptainer     - Build Apptainer image (alt to singularity)"
	@echo "  make run-sif       - Run SIF interactively"
	@echo "  make clean         - Remove Docker image & SIF file"

# ==== DOCKER ====
docker:
	@echo "🔧 Building Docker image: $(DOCKER_IMAGE)"
	docker build -t $(DOCKER_IMAGE) .

run:
	@echo "🚀 Running $(APP_NAME) on http://localhost:$(PORT)"
	docker run --rm \
		--name $(APP_NAME) \
		-p $(PORT):$(SHINY_PORT) \
		$(DOCKER_IMAGE)

run-dev:
	@echo "🧪 Running with local bind-mount (dev mode)"
	docker run --rm \
		--name $(APP_NAME)-dev \
		-p $(PORT):$(SHINY_PORT) \
		-v $(PWD):/app \
		$(DOCKER_IMAGE)

# ==== SINGULARITY ====
sif:
	@echo "📦 Building Singularity image: $(SIF_NAME)"
	singularity build $(SIF_NAME) Singularity.def

apptainer:
	@echo "📦 Building Apptainer image: $(SIF_NAME)"
	apptainer build $(SIF_NAME) Singularity.def

run-sif:
	@echo "🔥 Running Singularity image interactively"
	singularity run --bind $(PWD):/mnt $(SIF_NAME)

# ==== CLEANUP ====
clean:
	@echo "🧹 Removing Docker image and SIF"
	docker rmi $(DOCKER_IMAGE) || true
	rm -f $(SIF_NAME)
