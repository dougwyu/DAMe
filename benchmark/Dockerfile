FROM rust:1.85-bookworm@sha256:e51d0265072d2d9d5d320f6a44dde6b9ef13653b035098febd68cce8fa7c0bc4

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        python3=3.11.2-1+b1 \
        python3-pip=23.0.1+dfsg-1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY benchmark/constraints.txt /app/benchmark/constraints.txt
RUN python3 -m pip install --break-system-packages \
    --constraint /app/benchmark/constraints.txt numpy

COPY . /app

RUN cargo build --locked --release --manifest-path /app/rust/Cargo.toml
RUN python3 -m pip install --break-system-packages --no-deps /app/python

ENTRYPOINT ["python3", "/app/benchmark/run_v3_benchmarks.py"]
