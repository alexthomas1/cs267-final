# cs267-final
# Alex Thomas, Seth Park

In order to run unoptimized code input:
cargo run

To run the optimized code input:
cargo run --release

The different branches have our different implementation: <br />
master - Serial N-Body in Rust <br />
rayon - Rayon N-Body Implementation  <br />
crossbeam - Crossbeam N-Body Implementation  <br />
std_thread - Local buffer with single lock using the Rust std library  <br />
std_thread_lock - Every particle has a lock using the Rust std library  <br />
