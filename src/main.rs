use std::num;

extern crate rand;
extern crate clap;

use rand::Rng;
use std::{mem, ptr};
use std::time::{Duration, Instant};

const density: f64 = 0.0005;
const mass: f64 = 0.01;
const cutoff: f64 = 0.01;
const min_r: f64 = (cutoff) / 100.0;
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
fn set_size(n: i32)
{
    unsafe {
        size = (density * n as f64).sqrt();
    }
}

//
//  Initialize the particle positions and velocities
//


/* Bin each particle by pos
 *
 *
 * */
fn init_particles(n: i32, p: &mut [particle_t])
{
//    srand48(time(NULL));

    let mut rng = rand::thread_rng();

    let sx: i32 = ((n as f64).sqrt()).ceil() as i32;
    let sy: i32 = (n + sx - 1) / sx;

    let mut i: usize;
    let mut j: u32;
    let mut k: i32;

    let mut shuffle = vec![0; n as usize];

    for i in 0..n {
        shuffle[i as usize] = i;
    }

    for i in 0..n {

        //
        //  make sure particles are not spatially sorted
        //

        j = (rng.gen::<u32>() % ((n - i) as u32) ) as u32; //Generate random number from 0.. n-1
        k = shuffle[j as usize];
        shuffle[j as usize] = shuffle[(n - i - 1) as usize];
        //
        //  distribute particles evenly to ensure proper spacing
        //

        unsafe {
            p[i as usize].x = size * (1. + (k % sx) as f64) / (1.0 + sx as f64); //gets x
            p[i as usize].y = size * (1. + (k / sx) as f64) / (1.0 + sy as f64); //gets y
        }

        //
        //  assign random velocities within a bound
        //
        p[i as usize].vx = (rng.gen::<f64>() * 2.0 - 1.0) as f64;
        p[i as usize].vy = (rng.gen::<f64>() * 2.0 - 1.0) as f64;
    }
}

fn apply_force(p: &mut [particle_t], i: i32, j: i32, dmin: &mut f64, davg: &mut f64, navg: &mut i32)
{

    let dx: f64 = p[j as usize].x - p[i as usize].x;
    let dy: f64 = p[j as usize].y - p[i as usize].y;
    let mut r2: f64 = dx * dx + dy * dy;
    if r2 > cutoff * cutoff {
        return;
    }
    if r2 != 0.0
    {
        if r2 / (cutoff * cutoff) < *dmin * (*dmin) {
            *dmin = r2.sqrt() / cutoff;
        }
        (*davg) += r2.sqrt() / cutoff;
        (*navg) += 1;
    }

    r2 = r2.max(min_r * min_r);
    let r: f64 = r2.sqrt();


    //
    //  very simple short-range repulsive force
    //
    let coef: f64 = (1.0 - cutoff / r) / r2 / mass;
    p[i as usize].ax += coef * dx;
    p[i as usize].ay += coef * dy;
}


//
//  integrate the ODE
//

fn move_particle(p: &mut particle_t)
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    //
    //  bounce from walls
    //
    unsafe {
        while p.x < 0.0 || p.x > size {
            if p.x < 0.0 {
                p.x = -1.0 * p.x;
            } else {
                p.x = 2.0 * size - p.x;
            }
            p.vx = -p.vx;
        }
        while p.y < 0.0 || p.y > size {
            if p.y < 0.0 {
                p.y = -1.0 * p.y;
            } else {
                p.y = 2.0 * size - p.y;
            }
            p.vy = -p.vy;
        }
    }
}


fn main() {
    let mut navg: i32 = 0;
    let mut nabsavg: i32 = 0;
    const n: i32 = 1000;

    let mut davg: f64 = 0.0;
    let mut dmin: f64 = 0.0;
    let mut absmin: f64 = 1.0;
    let mut absavg: f64 = 0.0;


    let mut particles: [particle_t; n as usize] = unsafe { ::std::mem::uninitialized() };


    set_size(n);
    init_particles(n, &mut particles);


    let sx: i32 = ((n as f64).sqrt()).ceil() as i32;
    let sy: i32 = (n + sx - 1) / sx;

    let old_x: f64 = 0.0;
    let old_y: f64 = 0.0;
    let new_x: f64 = 0.0;
    let new_y: f64 = 0.0;

    let BOX_NUM: i32;
    unsafe { BOX_NUM = (size / cutoff).floor() as i32; }

    let mut step: i32 = 0;
    let mut i: usize = 0;
    let mut j: usize = 0;

    let now = Instant::now();

    for step in 0..NSTEPS {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;

        //
        //  compute forces
        //

        for i in 0..n {
            particles[i as usize].ax = 0.0;
            particles[i as usize].ay = 0.0;

            for j in 0..n {
                apply_force(&mut particles, i, j, &mut dmin, &mut davg, &mut navg);
            }
        }

        //
        //  move particles
        //
        for i in 0..n {
            move_particle (&mut particles[i as usize]);
//            println!("x: {}, y: {}", particles[i as usize].x, particles[i as usize].y);
        }


        //
        // Computing statistical data
        //
        if navg > 0 {
            absavg +=  davg/(navg as f64);
            nabsavg+=1;
        }

        if dmin < absmin {
            absmin = dmin;
        }


    }

    println!("n = {}, simulation time = {} seconds", n, now.elapsed().as_secs());

    if nabsavg > 0 {
        absavg /= nabsavg as f64;
    }

    //
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //

    println!(", absmin = {}, absavg = {}", absmin, absavg);
    if absmin < 0.4 {
        println!("The minimum distance is below 0.4 meaning that some particle is not interacting");
    }
    if absavg < 0.8 {
        println!("The average distance is below 0.8 meaning that most particles are not interacting");
    }

}






