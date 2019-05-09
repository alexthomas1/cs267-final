# cs267-final
# Alex Thomas, Seth Park

In order to run unoptimized code input:
cargo run

To run the optimized code input:
cargo run --release

The different branches have our different implementation:
master - Serial N-Body in Rust
rayon - Rayon N-Body Implementation
crossbeam - Crossbeam N-Body Implementation
std_thread - Local buffer with single lock using the Rust std library
std_thread_lock - Every particle has a lock using the Rust std library
