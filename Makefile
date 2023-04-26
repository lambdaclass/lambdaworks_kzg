.PHONY: test clippy

test:
	cargo test

clippy:
	cargo clippy --all-targets -- -D warnings
