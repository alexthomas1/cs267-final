const density: f64 = 0.0005;
const mass: f64 =  0.01;
const cutoff:f64 =0.01;
const min_r: f64 =  (cutoff/100);
const dt: f64 = 0.0005;

const NSTEPS: i32 = 1000;
const SAVEFREQ: i32 = 10;

static mut size: f64 = 0.0;


struct particle_t
{
     x: f64,
     y: f64,
     vx: f64,
     vy: f64,
     ax: f64,
     ay: f64,
     bin: f64,
}

//
//  keep density constant
//
fn set_size( n:i32 )
{
    size = (density * n).sqrt();
}

//
//  Initialize the particle positions and velocities
//


/* Bin each particle by pos
 *
 *
 * */
fn init_particles(n: i32, mut p: &particle_t)
{
    srand48(time(NULL));

    sx: f64 = ((n as f64).sqrt()).ceil();
    sy: f64  = (n as f64 + sx - 1)/sx;

    let mut shuffle:[i32; n] = [0; n];

    for i in 0..n{
        shuffle[i as usize] = i;
    }

    for i in 0..n{

        //
        //  make sure particles are not spatially sorted
        //
        j: i32 = lrand48() % (n - i);
        k: i32 = shuffle[j as usize];
        shuffle[j] = shuffle[n as usize - i as usize - 1];

        //
        //  distribute particles evenly to ensure proper spacing
        //

        p[i as usize].0 = size * (1. + (k % sx)) / (1 + sx); //gets x
        p[i as usize].1 = size * (1. + (k / sx)) / (1 + sy); //gets y

        //
        //  assign random velocities within a bound
        //
                p[i].2 = drand48() * 2 - 1;
                p[i].3 = drand48() * 2 - 1;
    }
}