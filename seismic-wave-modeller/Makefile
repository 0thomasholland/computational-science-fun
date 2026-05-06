CONFIG   ?= minimal_config.toml
ITER     ?= 0
FIELD    ?= vmag
FPS      ?= 30
BIN      := ./target/release/seismic-wave-modeller

.PHONY: build run run-example clean clean-output video gif help

## Build optimised binary
build:
	cargo build --release

## Run simulation (CONFIG=minimal_config.toml ITER=0)
run: build
	$(BIN) --config $(CONFIG) --iteration $(ITER) --verbose

## Run the full example config (3 sources, 5 km domain)
run-example: build
	$(BIN) --config example_config.toml --iteration 1 --verbose

## Stitch frames into mp4 then gif (ITER=0 FPS=30)
video:
	ffmpeg -y -framerate $(FPS) \
	       -pattern_type glob -i "output/$(ITER)/$(FIELD)_*.png" \
	       -c:v libx264 -pix_fmt yuv420p \
	       output/$(ITER).mp4
	ffmpeg -y -i output/$(ITER).mp4 \
	       -vf "fps=15,scale=720:-1:flags=lanczos" \
	       -c:v gif output/$(ITER).gif
	@echo "Written: output/$(ITER).mp4  output/$(ITER).gif"

## Remove build artefacts
clean:
	cargo clean

## Remove simulation output frames and videos
clean-output:
	rm -rf output/

## Show available targets
help:
	@grep -E '^##' Makefile | sed 's/## //'
